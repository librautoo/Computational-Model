%% Automatic three-delay atom splitting for every successful trial in one channel
% MATLAB R2024a; requires Optimization Toolbox (fmincon and lsqnonneg).
%
% Input
%   One original single-trial model-fit MAT-file for one EEG channel.
%
% Main output
%   1) summary.csv
%   2) noise_stability_summary.csv
%
% Scientific role of this script
%   The number of atoms is fixed at K = 3. The script estimates only the
%   three continuous delay locations as the scientific output of this
%   stage. It does not compare K = 1, 2, and 3, and it does not use BIC to
%   choose the number of atoms.
%
%   For numerical evaluation of a candidate delay vector, the nonnegative
%   atom coefficients and a nonnegative baseline are solved internally by
%   least squares:
%
%       min_delta  min_{a >= 0, b >= 0}
%           ||M(t) - b - sum_v a_v S0(t-delta_v)||_2^2.
%
%   These coefficients are nuisance quantities only. They are not the
%   final physiological proportions and are not a replacement for the
%   later conditional re-estimation of theta. After the three delays are
%   accepted, fix them and estimate theta in the full three-delay model.
%
% Atom-splitting strategy
%   A one-atom solution is fitted first. The largest active atom is split
%   into two nearby atoms and the delays are jointly refined. The largest
%   atom is then split again, and the final three delays are jointly
%   refined in continuous time. Linear nuisance coefficients are profiled
%   out at every objective evaluation.
%
% Noise stability
%   The complete three-delay estimator is repeated after adding Gaussian
%   white noise with standard deviations equal to 1%, 2%, and 5% of the
%   trial-wise observed-signal standard deviation. Ten reproducible runs
%   are performed at each level. A noisy result is delay-stable when all
%   three sorted delays remain within the configured tolerance of the
%   unperturbed result.

clearvars;
clc;

%% USER SETTINGS
cfg = struct;
cfg.dataDir = fileparts(mfilename('fullpath'));

% TO RUN ANOTHER DATASET/CHANNEL, CHANGE ONLY THIS FILE NAME.
% Keep the exact suffix, including "(1)", when it is part of the filename.
cfg.fitFile = ...
    'DataA_Fz_single_trial_original_model_fit_fixed_ranges.mat';

% Results are written to:
% <input folder>/atom_splitting_results/Data_A-or-Data_B/<channel>/
cfg.outputFolder = 'atom_splitting_results';

% Both experimental directions are processed by default.
cfg.processDirections = ["Right", "Left"];

% The scientific model order is fixed. Do not change this value for the
% intended three-delay analysis.
cfg.numberOfAtoms = 3;

% Delay constraints. At 256 Hz, one sample is 3.90625 ms.
cfg.minimumDelaySeparationMs = 3.90625;
cfg.edgeMarginMs = 0;

% Each split places two child atoms on opposite sides of the parent before
% joint continuous-delay refinement. Three samples at 256 Hz = 11.71875 ms.
cfg.splitHalfWidthMs = 11.71875;

% Multistart settings for K = 1, 2, and 3. The deterministic split seeds are
% always included; these numbers add reproducible random starts.
cfg.randomStartsPerStage = [2, 4, 8];
cfg.randomJitterMs = 15.625; % four samples at 256 Hz
cfg.optimizationBaseSeed = 314159;

% This threshold is only an identifiability diagnostic. If an internally
% fitted nuisance coefficient contributes less than 1% of their total,
% that atom is effectively inactive and its delay is not considered
% identifiable. The shares are not interpreted as physiological ratios.
cfg.minimumAuxiliaryShare = 0.01;

% Keep false for consistency with the transition-centered S0 proxy used by
% the current MUSIC workflow. Setting true is a sensitivity analysis only.
cfg.removeTemplateFloor = false;

% Optimizer settings.
cfg.optimizerMaxIterations = 800;
cfg.optimizerMaxFunctionEvaluations = 12000;
cfg.optimizerFunctionTolerance = 1e-10;
cfg.optimizerStepTolerance = 1e-11;

% Noise-stability settings: 1%, 2%, and 5%, with 10 runs per level.
cfg.runNoiseStability = true;
cfg.noiseStabilityLevels = [1e-2, 2e-2, 5e-2];
cfg.noiseStabilityReplicates = 10;
cfg.noiseStabilityBaseSeed = 271828;
cfg.noiseStabilityDelayToleranceMs = 7.8125; % two samples at 256 Hz

%% Validate the environment and load one channel-level file
if exist('fmincon', 'file') ~= 2 || exist('lsqnonneg', 'file') ~= 2
    error(['This script requires fmincon and lsqnonneg from ' ...
        'Optimization Toolbox.']);
end

if cfg.numberOfAtoms ~= 3
    error('cfg.numberOfAtoms must remain equal to 3 for this analysis.');
end
if numel(cfg.randomStartsPerStage) ~= cfg.numberOfAtoms
    error('cfg.randomStartsPerStage must contain one value for each stage.');
end

inputPath = fullfile(cfg.dataDir, cfg.fitFile);
if ~isfile(inputPath)
    error(['Input fit file was not found:\n  %s\n' ...
        'Change only cfg.fitFile to the exact MAT-file name.'], inputPath);
end

loaded = load(inputPath, 'FitAnalysis');
if ~isfield(loaded, 'FitAnalysis')
    error('The selected MAT-file does not contain FitAnalysis.');
end
fitAnalysis = loaded.FitAnalysis;
if ~isfield(fitAnalysis, 'finalResults')
    error('FitAnalysis does not contain finalResults.');
end
if ~isfield(fitAnalysis, 'model') || ...
        ~isfield(fitAnalysis.model, 'parameter_names')
    error('FitAnalysis.model.parameter_names is missing.');
end

if isfield(fitAnalysis, 'dataset')
    datasetName = scalar_string(fitAnalysis.dataset);
else
    datasetName = "unknown_dataset";
end
if isfield(fitAnalysis, 'channel')
    channelName = scalar_string(fitAnalysis.channel);
else
    channelName = "unknown_channel";
end

parameterNames = string(fitAnalysis.model.parameter_names(:));
t1Index = find(strcmpi(parameterNames, 't1'), 1);
if isempty(t1Index)
    error('The original fitted parameter t1 was not found.');
end

allResults = fitAnalysis.finalResults;
nAll = numel(allResults);
directions = strings(nAll, 1);
trialNumbers = nan(nAll, 1);
successful = false(nAll, 1);
for resultIndex = 1:nAll
    directions(resultIndex) = ...
        scalar_string(allResults(resultIndex).direction);
    trialNumbers(resultIndex) = ...
        double(allResults(resultIndex).trial_number);
    successful(resultIndex) = ...
        logical(allResults(resultIndex).success);
end

% Process Right trials and then Left trials, each in trial-number order.
selectedIndices = [];
for requestedDirection = cfg.processDirections
    candidates = find(successful & ...
        strcmpi(directions, requestedDirection));
    [~, chronologicalOrder] = sort(trialNumbers(candidates));
    selectedIndices = [selectedIndices; ...
        candidates(chronologicalOrder)]; %#ok<AGROW>
end
selectedIndices = unique(selectedIndices, 'stable');
if isempty(selectedIndices)
    error('No successful trials matched cfg.processDirections.');
end

fprintf('\nDataset: %s | Channel: %s\n', ...
    scalar_char(datasetName), scalar_char(channelName));
fprintf('Processing %d successful trials with fixed K = 3.\n', ...
    numel(selectedIndices));

%% Create the same DataA/DataB and channel folder hierarchy as MUSIC
safeChannel = string(regexprep(channelName, ...
    '[^A-Za-z0-9_-]', '_'));
[~, sourceStemForFolder] = fileparts(inputPath);
sourceStemForFolder = string(sourceStemForFolder);

if startsWith(sourceStemForFolder, "DataA_", 'IgnoreCase', true)
    outputDatasetFolder = "Data_A";
elseif startsWith(sourceStemForFolder, "DataB_", 'IgnoreCase', true)
    outputDatasetFolder = "Data_B";
else
    outputDatasetFolder = string(regexprep(datasetName, ...
        '[^A-Za-z0-9_-]', '_'));
    warning(['Could not determine DataA/DataB from the input filename. ' ...
        'Using FitAnalysis.dataset for the output folder instead.']);
end

outputDir = fullfile(fileparts(inputPath), cfg.outputFolder, ...
    outputDatasetFolder, safeChannel);
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

summaryRows = repmat(empty_summary_row(), numel(selectedIndices), 1);

%% Estimate three delays for every successful trial
for trialPosition = 1:numel(selectedIndices)
    result = allResults(selectedIndices(trialPosition));
    direction = scalar_string(result.direction);
    trialNumber = double(result.trial_number);
    optimizerSeed = trial_optimizer_seed( ...
        cfg, direction, trialNumber);

    fprintf('\n[%d/%d] %s trial %g\n', ...
        trialPosition, numel(selectedIndices), ...
        scalar_char(direction), trialNumber);

    summaryRows(trialPosition) = safe_estimate_trial( ...
        result, datasetName, channelName, t1Index, ...
        [], 0, NaN, optimizerSeed, cfg);

    row = summaryRows(trialPosition);
    if row.ThreeDelayValid
        fprintf(['  FOUND: delta = [%.3f, %.3f, %.3f] ms | ' ...
            'R^2 = %.4f | min auxiliary share = %.3f\n'], ...
            row.Delta1Ms, row.Delta2Ms, row.Delta3Ms, ...
            row.ProfiledR2, row.MinimumAuxiliaryShare);
    elseif row.OptimizationSucceeded
        fprintf(['  Candidate: delta = [%.3f, %.3f, %.3f] ms | ' ...
            'status = %s\n'], ...
            row.Delta1Ms, row.Delta2Ms, row.Delta3Ms, ...
            scalar_char(row.Status));
    else
        fprintf('  FAILED: %s\n', scalar_char(row.Message));
    end
end

%% Run the matched noise-stability experiment
noiseStabilitySummaryTable = struct2table( ...
    repmat(empty_stability_summary_row(), 0, 1));
if cfg.runNoiseStability
    fprintf('\nStarting noise stability: %d levels x %d replicates.\n', ...
        numel(cfg.noiseStabilityLevels), ...
        cfg.noiseStabilityReplicates);

    noiseStabilitySummaryTable = run_noise_stability( ...
        allResults, selectedIndices, summaryRows, ...
        datasetName, channelName, t1Index, cfg);
end

%% Save channel-level output files
summaryTable = struct2table(summaryRows);
[~, sourceFitStem, sourceFitExtension] = fileparts(inputPath);
sourceFitFile = string([sourceFitStem sourceFitExtension]);
summaryTable.SourceFitFile = repmat( ...
    sourceFitFile, height(summaryTable), 1);

summaryCsvPath = fullfile(outputDir, 'summary.csv');
stabilityCsvPath = fullfile(outputDir, ...
    'noise_stability_summary.csv');
writetable(summaryTable, summaryCsvPath);
if cfg.runNoiseStability
    noiseStabilitySummaryTable.SourceFitFile = repmat( ...
        sourceFitFile, height(noiseStabilitySummaryTable), 1);
    writetable(noiseStabilitySummaryTable, stabilityCsvPath);
end

validCount = nnz(summaryTable.ThreeDelayValid);
fprintf('\n============================================================\n');
fprintf('Channel run completed: %d/%d trials have identifiable three delays.\n', ...
    validCount, height(summaryTable));
fprintf('Delay summary:      %s\n', summaryCsvPath);
if cfg.runNoiseStability
    fprintf('Noise stability:    %s\n', stabilityCsvPath);
end
fprintf('============================================================\n');

displayColumns = {'Direction', 'TrialNumber', 'ThreeDelayValid', ...
    'Delta1Ms', 'Delta2Ms', 'Delta3Ms', ...
    'MinimumSeparationMs', 'ThreeAtomsActive', ...
    'ProfiledR2', 'Status'};
disp(summaryTable(:, displayColumns));

%% Local function: safely estimate one trial
function row = safe_estimate_trial( ...
        result, datasetName, channelName, t1Index, ...
        observationOverride, noiseStd, noiseSeed, optimizerSeed, cfg)
    try
        row = estimate_trial( ...
            result, datasetName, channelName, t1Index, ...
            observationOverride, noiseStd, noiseSeed, ...
            optimizerSeed, cfg);
    catch ME
        row = empty_summary_row();
        row.Dataset = datasetName;
        row.Channel = channelName;
        row.Direction = scalar_string(result.direction);
        row.TrialNumber = double(result.trial_number);
        row.NoiseStd = noiseStd;
        row.NoiseSeed = noiseSeed;
        row.OptimizerSeed = optimizerSeed;
        row.Status = "trial_error";
        row.Message = string(ME.message);
    end
end

%% Local function: estimate fixed K = 3 delays for one observation
function row = estimate_trial( ...
        result, datasetName, channelName, t1Index, ...
        observationOverride, noiseStd, noiseSeed, optimizerSeed, cfg)

    row = empty_summary_row();
    row.Dataset = datasetName;
    row.Channel = channelName;
    row.Direction = scalar_string(result.direction);
    row.TrialNumber = double(result.trial_number);
    row.NoiseStd = noiseStd;
    row.NoiseSeed = noiseSeed;
    row.OptimizerSeed = optimizerSeed;

    if isfield(result, 'N_coup')
        row.OriginalNcoup = double(result.N_coup);
    end
    if isfield(result, 'metrics') && isfield(result.metrics, 'R2')
        row.OriginalR2 = double(result.metrics.R2);
    end

    tAll = double(result.t_data(:));
    observedOriginalAll = double(result.y_data(:));
    fittedAll = double(result.y_model(:));
    parameters = double(result.parameters(:));
    if numel(parameters) < t1Index
        error('The result parameter vector does not contain t1.');
    end
    oldT1 = parameters(t1Index);

    if isempty(observationOverride)
        observedAll = observedOriginalAll;
    else
        observedAll = double(observationOverride(:));
        if numel(observedAll) ~= numel(tAll)
            error('The noisy observation length does not match t_data.');
        end
    end

    validMask = isfinite(tAll) & isfinite(observedAll) & ...
        isfinite(observedOriginalAll) & isfinite(fittedAll);
    t = tAll(validMask);
    observed = observedAll(validMask);
    fitted = fittedAll(validMask);

    if numel(t) < 8
        error('Fewer than eight finite samples remain in this trial.');
    end
    if any(diff(t) <= 0)
        error('t_data must be strictly increasing.');
    end
    if std(observed) <= eps
        error('The observed response has no usable variation.');
    end

    dt = median(diff(t));
    fitStart = t(1) + cfg.edgeMarginMs/1000;
    fitEnd = t(end) - cfg.edgeMarginMs/1000;
    minimumSeparation = cfg.minimumDelaySeparationMs/1000;
    if fitEnd-fitStart < ...
            (cfg.numberOfAtoms-1)*minimumSeparation
        error('The fit window is too short for three separated delays.');
    end

    localTime = t-oldT1;
    if min(localTime) >= 0 || max(localTime) <= 0
        error('The original fitted t1 is outside the trial fit window.');
    end

    atomTemplate = fitted;
    if cfg.removeTemplateFloor
        atomTemplate = atomTemplate-min(atomTemplate);
    end
    templateScale = max(abs(atomTemplate));
    if ~isfinite(templateScale) || templateScale <= eps
        error('The fitted base-response proxy has invalid amplitude.');
    end
    atomTemplate = atomTemplate/templateScale;

    row.OriginalT1Ms = 1000*oldT1;
    row.FsHz = 1/dt;
    row.FitStartMs = 1000*fitStart;
    row.FitEndMs = 1000*fitEnd;
    row.NoiseStdFraction = safe_noise_fraction( ...
        noiseStd, observedOriginalAll);

    stream = RandStream('mt19937ar', 'Seed', optimizerSeed);
    model = split_to_three_atoms( ...
        t, observed, localTime, atomTemplate, oldT1, ...
        fitStart, fitEnd, minimumSeparation, stream, cfg);

    row.OptimizationSucceeded = model.optimizationSucceeded;
    row.Converged = model.converged;
    row.ExitFlag = model.exitFlag;
    row.TotalStarts = model.totalAllStageStarts;
    row.SuccessfulStarts = model.successfulAllStageStarts;
    row.Delta1Ms = model.delaysMs(1);
    row.Delta2Ms = model.delaysMs(2);
    row.Delta3Ms = model.delaysMs(3);
    row.MinimumSeparationMs = min(diff(model.delaysMs));
    row.AllInsideFitWindow = all(model.delays >= fitStart-1e-10 & ...
        model.delays <= fitEnd+1e-10);
    row.WindowValid = row.AllInsideFitWindow && ...
        row.MinimumSeparationMs >= ...
        cfg.minimumDelaySeparationMs-1e-7;

    row.AuxiliaryCoefficient1 = model.amplitudes(1);
    row.AuxiliaryCoefficient2 = model.amplitudes(2);
    row.AuxiliaryCoefficient3 = model.amplitudes(3);
    row.MinimumAuxiliaryShare = min(model.auxiliaryShares);
    row.ThreeAtomsActive = all(model.auxiliaryShares >= ...
        cfg.minimumAuxiliaryShare);

    row.Baseline = model.baseline;
    row.ProfiledR2 = model.R2;
    row.ProfiledRMSE = model.RMSE;
    row.ProfiledSSE = model.SSE;

    % Positive convergence is reported separately. A finite optimized
    % solution can still be inspected when fmincon stops at its iteration
    % limit, but it is not automatically labeled as a valid final result.
    row.ThreeDelayValid = row.OptimizationSucceeded && ...
        row.Converged && row.WindowValid && row.ThreeAtomsActive;

    if row.ThreeDelayValid
        row.Status = "ok_three_delays";
    elseif ~row.Converged
        row.Status = "optimizer_not_converged";
    elseif ~row.WindowValid
        row.Status = "delay_constraints_not_satisfied";
    elseif ~row.ThreeAtomsActive
        row.Status = "inactive_auxiliary_atom";
    else
        row.Status = "invalid_three_delay_result";
    end
end

%% Local function: grow K = 1 into the fixed K = 3 model
function finalModel = split_to_three_atoms( ...
        t, y, localTime, atomTemplate, oldT1, ...
        fitStart, fitEnd, minimumSeparation, stream, cfg)

    models = repmat(empty_model_result(), cfg.numberOfAtoms, 1);

    for numberOfAtoms = 1:cfg.numberOfAtoms
        lowerBounds = fitStart*ones(1, numberOfAtoms);
        upperBounds = fitEnd*ones(1, numberOfAtoms);

        % Include several deterministic seeds. The first seed follows the
        % actual splitting path; the others protect against a poor parent
        % atom without changing the fixed K = 3 scientific target.
        if numberOfAtoms == 1
            [~, peakIndex] = max(y);
            deterministicSeeds = [ ...
                oldT1; ...
                t(peakIndex); ...
                (fitStart+fitEnd)/2];
        else
            splitSeed = split_largest_atom( ...
                models(numberOfAtoms-1), ...
                cfg.splitHalfWidthMs/1000, ...
                lowerBounds, upperBounds, minimumSeparation);

            equallySpaced = linspace( ...
                fitStart, fitEnd, numberOfAtoms+2);
            equallySpaced = equallySpaced(2:end-1);

            centered = oldT1 + ...
                ((1:numberOfAtoms)-(numberOfAtoms+1)/2) * ...
                (2*cfg.splitHalfWidthMs/1000);

            deterministicSeeds = [ ...
                splitSeed; equallySpaced; centered];
        end

        for seedIndex = 1:size(deterministicSeeds, 1)
            deterministicSeeds(seedIndex, :) = make_feasible_delays( ...
                deterministicSeeds(seedIndex, :), ...
                lowerBounds, upperBounds, minimumSeparation);
        end

        models(numberOfAtoms) = fit_profiled_atomic_model( ...
            t, y, localTime, atomTemplate, numberOfAtoms, ...
            deterministicSeeds, lowerBounds, upperBounds, ...
            minimumSeparation, ...
            cfg.randomStartsPerStage(numberOfAtoms), stream, cfg);
    end

    finalModel = models(end);
    finalModel.totalAllStageStarts = sum([models.totalStarts]);
    finalModel.successfulAllStageStarts = ...
        sum([models.successfulStarts]);
end

%% Local function: fit continuous delays while profiling out coefficients
function model = fit_profiled_atomic_model( ...
        t, y, localTime, atomTemplate, numberOfAtoms, ...
        deterministicSeeds, lowerBounds, upperBounds, ...
        minimumSeparation, randomStartCount, stream, cfg)

    deterministicSeeds = double(deterministicSeeds);
    if isvector(deterministicSeeds)
        deterministicSeeds = reshape( ...
            deterministicSeeds, [], numberOfAtoms);
    end
    if size(deterministicSeeds, 2) ~= numberOfAtoms
        error('A deterministic delay seed has the wrong width.');
    end

    lowerBounds = double(lowerBounds(:)).';
    upperBounds = double(upperBounds(:)).';

    Aineq = zeros(max(numberOfAtoms-1, 0), numberOfAtoms);
    bineq = -minimumSeparation * ...
        ones(max(numberOfAtoms-1, 0), 1);
    for delayIndex = 1:numberOfAtoms-1
        Aineq(delayIndex, delayIndex) = 1;
        Aineq(delayIndex, delayIndex+1) = -1;
    end

    options = optimoptions('fmincon', ...
        'Display', 'off', ...
        'Algorithm', 'sqp', ...
        'MaxIterations', cfg.optimizerMaxIterations, ...
        'MaxFunctionEvaluations', ...
        cfg.optimizerMaxFunctionEvaluations, ...
        'FunctionTolerance', cfg.optimizerFunctionTolerance, ...
        'OptimalityTolerance', 1e-9, ...
        'StepTolerance', cfg.optimizerStepTolerance, ...
        'FiniteDifferenceType', 'central');

    deterministicCount = size(deterministicSeeds, 1);
    totalStarts = deterministicCount + randomStartCount;
    bestObjective = inf;
    bestDelays = [];
    bestExitFlag = NaN;
    bestOutput = struct;
    successfulStarts = 0;

    objective = @(delays) profiled_atomic_sse( ...
        delays, t, y, localTime, atomTemplate);

    for startIndex = 1:totalStarts
        if startIndex <= deterministicCount
            initialDelays = deterministicSeeds(startIndex, :);
        elseif mod(startIndex-deterministicCount, 2) == 1
            referenceIndex = 1+mod( ...
                startIndex-deterministicCount-1, deterministicCount);
            initialDelays = deterministicSeeds(referenceIndex, :) + ...
                (cfg.randomJitterMs/1000) * ...
                randn(stream, 1, numberOfAtoms);
        else
            initialDelays = lowerBounds + ...
                rand(stream, 1, numberOfAtoms) .* ...
                (upperBounds-lowerBounds);
        end

        initialDelays = make_feasible_delays( ...
            initialDelays, lowerBounds, upperBounds, ...
            minimumSeparation);

        try
            [fittedDelays, objectiveValue, exitFlag, output] = ...
                fmincon(objective, initialDelays, ...
                Aineq, bineq, [], [], lowerBounds, upperBounds, ...
                [], options);
        catch
            continue;
        end

        if isfinite(objectiveValue)
            successfulStarts = successfulStarts+1;
        end
        if isfinite(objectiveValue) && objectiveValue < bestObjective
            bestObjective = objectiveValue;
            bestDelays = fittedDelays;
            bestExitFlag = exitFlag;
            bestOutput = output;
        end
    end

    if isempty(bestDelays)
        error('All optimization starts failed for K = %d.', numberOfAtoms);
    end

    delays = sort(double(bestDelays(:)).');
    [amplitudes, baseline, atomMatrix, fittedSignal] = ...
        profiled_linear_fit(t, y, localTime, atomTemplate, delays);

    residual = y(:)-fittedSignal(:);
    SSE = sum(residual.^2);
    SST = sum((y(:)-mean(y(:))).^2);
    if SST > eps
        R2 = 1-SSE/SST;
    else
        R2 = NaN;
    end

    amplitudeSum = sum(amplitudes);
    if amplitudeSum > eps
        auxiliaryShares = amplitudes/amplitudeSum;
    else
        auxiliaryShares = zeros(size(amplitudes));
    end

    model = empty_model_result();
    model.numberOfAtoms = numberOfAtoms;
    model.delays = delays;
    model.delaysMs = 1000*delays;
    model.amplitudes = amplitudes;
    model.auxiliaryShares = auxiliaryShares;
    model.baseline = baseline;
    model.atomMatrix = atomMatrix;
    model.fittedSignal = fittedSignal;
    model.residual = residual;
    model.SSE = SSE;
    model.R2 = R2;
    model.RMSE = sqrt(mean(residual.^2));
    model.exitFlag = bestExitFlag;
    model.converged = bestExitFlag > 0;
    model.optimizationSucceeded = isfinite(SSE);
    model.iterations = get_output_value( ...
        bestOutput, 'iterations', NaN);
    model.functionEvaluations = get_output_value( ...
        bestOutput, 'funcCount', NaN);
    model.totalStarts = totalStarts;
    model.successfulStarts = successfulStarts;
end

%% Local function: variable-projection objective
function SSE = profiled_atomic_sse( ...
        delays, t, y, localTime, atomTemplate)
    try
        [~, ~, ~, fittedSignal] = profiled_linear_fit( ...
            t, y, localTime, atomTemplate, delays);
        residual = y(:)-fittedSignal(:);
        if any(~isfinite(residual))
            SSE = realmax('double');
        else
            SSE = sum(residual.^2);
        end
    catch
        SSE = realmax('double');
    end
end

%% Local function: solve nuisance coefficients for fixed delays
function [amplitudes, baseline, atomMatrix, fittedSignal] = ...
        profiled_linear_fit(t, y, localTime, atomTemplate, delays)
    atomMatrix = build_atom_matrix( ...
        t, localTime, atomTemplate, delays);
    designMatrix = [atomMatrix, ones(numel(t), 1)];
    coefficients = lsqnonneg(designMatrix, y(:));
    amplitudes = coefficients(1:end-1).';
    baseline = coefficients(end);
    fittedSignal = baseline + atomMatrix*amplitudes(:);
end

%% Local function: construct continuously shifted atom columns
function atomMatrix = build_atom_matrix( ...
        t, localTime, atomTemplate, delays)
    t = double(t(:));
    localTime = double(localTime(:));
    atomTemplate = double(atomTemplate(:));
    delays = double(delays(:)).';

    atomMatrix = zeros(numel(t), numel(delays));
    for delayIndex = 1:numel(delays)
        queryLocalTime = t-delays(delayIndex);

        % Zero extrapolation matches the zero-padded, transition-centered
        % S0 construction used by the current MUSIC script.
        atomMatrix(:, delayIndex) = interp1( ...
            localTime, atomTemplate, queryLocalTime, 'pchip', 0);
    end
    atomMatrix(~isfinite(atomMatrix)) = 0;
end

%% Local function: split the largest current atom
function splitDelays = split_largest_atom( ...
        previousModel, splitHalfWidth, lowerBounds, ...
        upperBounds, minimumSeparation)
    [~, largestIndex] = max(previousModel.amplitudes);
    oldDelays = previousModel.delays;
    splitCenter = oldDelays(largestIndex);
    oldDelays(largestIndex) = [];

    splitDelays = sort([oldDelays, ...
        splitCenter-splitHalfWidth, ...
        splitCenter+splitHalfWidth]);
    splitDelays = make_feasible_delays( ...
        splitDelays, lowerBounds, upperBounds, ...
        minimumSeparation);
end

%% Local function: enforce delay order and minimum separation
function delays = make_feasible_delays( ...
        delays, lowerBounds, upperBounds, minimumSeparation)
    delays = sort(double(delays(:)).');
    lowerBounds = double(lowerBounds(:)).';
    upperBounds = double(upperBounds(:)).';
    numberOfDelays = numel(delays);

    delays = min(max(delays, lowerBounds), upperBounds);
    for passIndex = 1:6 %#ok<NASGU>
        for delayIndex = 2:numberOfDelays
            delays(delayIndex) = max(delays(delayIndex), ...
                delays(delayIndex-1)+minimumSeparation);
        end
        delays = min(delays, upperBounds);

        for delayIndex = numberOfDelays-1:-1:1
            delays(delayIndex) = min(delays(delayIndex), ...
                delays(delayIndex+1)-minimumSeparation);
        end
        delays = max(delays, lowerBounds);
    end

    if any(delays < lowerBounds-1e-12) || ...
            any(delays > upperBounds+1e-12) || ...
            any(diff(delays) < minimumSeparation-1e-12)
        error('The requested delay bounds contain no feasible solution.');
    end
end

%% Local function: matched noise-stability experiment
function stabilityTable = run_noise_stability( ...
        allResults, selectedIndices, baselineRows, ...
        datasetName, channelName, t1Index, cfg)

    noiseLevels = double(cfg.noiseStabilityLevels(:)).';
    replicateCount = double(cfg.noiseStabilityReplicates);
    maximumRows = numel(selectedIndices)*numel(noiseLevels);
    outputRows = repmat(empty_stability_summary_row(), ...
        maximumRows, 1);
    outputCount = 0;

    for trialPosition = 1:numel(selectedIndices)
        result = allResults(selectedIndices(trialPosition));
        baseline = baselineRows(trialPosition);
        direction = scalar_string(result.direction);
        trialNumber = double(result.trial_number);
        optimizerSeed = trial_optimizer_seed( ...
            cfg, direction, trialNumber);

        fprintf('\n  Stability [%d/%d] %s trial %g\n', ...
            trialPosition, numel(selectedIndices), ...
            scalar_char(direction), trialNumber);

        for levelIndex = 1:numel(noiseLevels)
            noiseLevel = noiseLevels(levelIndex);
            runRows = repmat(empty_stability_run_row(), ...
                replicateCount, 1);

            for replicateIndex = 1:replicateCount
                noiseSeed = stability_noise_seed( ...
                    cfg, direction, trialNumber, ...
                    levelIndex, replicateIndex);
                [noisyObservation, noiseStd] = ...
                    make_noisy_observation( ...
                    result.y_data, noiseLevel, noiseSeed);

                noisy = safe_estimate_trial( ...
                    result, datasetName, channelName, t1Index, ...
                    noisyObservation, noiseStd, noiseSeed, ...
                    optimizerSeed, cfg);

                runRows(replicateIndex) = make_stability_run_row( ...
                    noisy, baseline, noiseLevel, replicateIndex, cfg);
            end

            outputCount = outputCount+1;
            outputRows(outputCount) = summarize_stability_runs( ...
                runRows, baseline, cfg);

            current = outputRows(outputCount);
            fprintf(['    noise=%5.2f%%: valid %.2f | ' ...
                'stable among valid %.2f | stable overall %.2f\n'], ...
                current.NoiseLevelPercent, ...
                current.ThreeDelayValidProportion, ...
                current.DelayStableProportionAmongValid, ...
                current.DelayStableProportionAllReplicates);
        end
    end

    stabilityTable = struct2table(outputRows(1:outputCount));
end

%% Local function: add reproducible Gaussian white noise
function [noisyObservation, noiseStd] = ...
        make_noisy_observation(observedOriginal, noiseLevel, noiseSeed)
    observedOriginal = double(observedOriginal(:));
    finiteValues = observedOriginal(isfinite(observedOriginal));
    if numel(finiteValues) < 2
        error('The original observation has too few finite samples.');
    end

    noiseStd = double(noiseLevel)*std(finiteValues);
    stream = RandStream('mt19937ar', 'Seed', noiseSeed);
    noise = noiseStd*randn(stream, size(observedOriginal));
    noisyObservation = observedOriginal+noise;
end

%% Local function: convert one noisy estimate to a stability row
function row = make_stability_run_row( ...
        noisy, baseline, noiseLevel, replicateIndex, cfg)
    row = empty_stability_run_row();
    row.Dataset = noisy.Dataset;
    row.Channel = noisy.Channel;
    row.Direction = noisy.Direction;
    row.TrialNumber = noisy.TrialNumber;
    row.NoiseLevelFraction = noiseLevel;
    row.NoiseLevelPercent = 100*noiseLevel;
    row.Replicate = replicateIndex;
    row.NoiseStd = noisy.NoiseStd;
    row.NoiseSeed = noisy.NoiseSeed;
    row.OptimizationSucceeded = noisy.OptimizationSucceeded;
    row.Converged = noisy.Converged;
    row.ThreeAtomsActive = noisy.ThreeAtomsActive;
    row.ThreeDelayValid = noisy.ThreeDelayValid;
    row.Delta1Ms = noisy.Delta1Ms;
    row.Delta2Ms = noisy.Delta2Ms;
    row.Delta3Ms = noisy.Delta3Ms;
    row.Status = noisy.Status;
    row.Message = noisy.Message;

    row.BaselineThreeDelayValid = baseline.ThreeDelayValid;
    if baseline.ThreeDelayValid && noisy.ThreeDelayValid
        baselineDelays = [baseline.Delta1Ms, ...
            baseline.Delta2Ms, baseline.Delta3Ms];
        noisyDelays = [noisy.Delta1Ms, ...
            noisy.Delta2Ms, noisy.Delta3Ms];
        changes = abs(noisyDelays-baselineDelays);
        row.Delta1AbsoluteChangeMs = changes(1);
        row.Delta2AbsoluteChangeMs = changes(2);
        row.Delta3AbsoluteChangeMs = changes(3);
        row.MaximumAbsoluteDelayChangeMs = max(changes);
        row.DelaysWithinToleranceOfBaseline = ...
            all(changes <= cfg.noiseStabilityDelayToleranceMs);
    end
end

%% Local function: summarize one trial and one noise level
function row = summarize_stability_runs(runRows, baseline, cfg)
    row = empty_stability_summary_row();
    row.Dataset = runRows(1).Dataset;
    row.Channel = runRows(1).Channel;
    row.Direction = runRows(1).Direction;
    row.TrialNumber = runRows(1).TrialNumber;
    row.NoiseLevelFraction = runRows(1).NoiseLevelFraction;
    row.NoiseLevelPercent = runRows(1).NoiseLevelPercent;
    row.ReplicateCount = numel(runRows);
    row.DelayToleranceMs = cfg.noiseStabilityDelayToleranceMs;
    row.BaselineThreeDelayValid = baseline.ThreeDelayValid;
    row.BaselineDelta1Ms = baseline.Delta1Ms;
    row.BaselineDelta2Ms = baseline.Delta2Ms;
    row.BaselineDelta3Ms = baseline.Delta3Ms;

    succeeded = [runRows.OptimizationSucceeded];
    converged = [runRows.Converged];
    active = [runRows.ThreeAtomsActive];
    valid = [runRows.ThreeDelayValid];
    stable = [runRows.DelaysWithinToleranceOfBaseline];

    row.OptimizationSucceededCount = nnz(succeeded);
    row.OptimizationSucceededProportion = mean(double(succeeded));
    row.ConvergedCount = nnz(converged);
    row.ConvergedProportion = mean(double(converged));
    row.ThreeAtomsActiveCount = nnz(active);
    row.ThreeAtomsActiveProportion = mean(double(active));
    row.ThreeDelayValidCount = nnz(valid);
    row.ThreeDelayValidProportion = mean(double(valid));
    row.DelayStableCount = nnz(stable);
    row.DelayStableProportionAllReplicates = ...
        row.DelayStableCount/row.ReplicateCount;
    if row.ThreeDelayValidCount > 0 && baseline.ThreeDelayValid
        row.DelayStableProportionAmongValid = ...
            row.DelayStableCount/row.ThreeDelayValidCount;
    end

    [row.Delta1MeanMs, row.Delta1MedianMs, ...
        row.Delta1StdMs, row.Delta1CI95LowMs, ...
        row.Delta1CI95HighMs] = finite_statistics( ...
        [runRows(valid).Delta1Ms]);
    [row.Delta2MeanMs, row.Delta2MedianMs, ...
        row.Delta2StdMs, row.Delta2CI95LowMs, ...
        row.Delta2CI95HighMs] = finite_statistics( ...
        [runRows(valid).Delta2Ms]);
    [row.Delta3MeanMs, row.Delta3MedianMs, ...
        row.Delta3StdMs, row.Delta3CI95LowMs, ...
        row.Delta3CI95HighMs] = finite_statistics( ...
        [runRows(valid).Delta3Ms]);

    comparable = valid & baseline.ThreeDelayValid;
    [row.MeanMaximumAbsoluteDelayChangeMs, ...
        row.MedianMaximumAbsoluteDelayChangeMs, ~, ~, ~] = ...
        finite_statistics( ...
        [runRows(comparable).MaximumAbsoluteDelayChangeMs]);
end

%% Local utility: finite descriptive statistics and empirical 95% interval
function [meanValue, medianValue, stdValue, ciLow, ciHigh] = ...
        finite_statistics(values)
    values = double(values(:));
    values = values(isfinite(values));
    if isempty(values)
        meanValue = NaN;
        medianValue = NaN;
        stdValue = NaN;
        ciLow = NaN;
        ciHigh = NaN;
        return;
    end

    meanValue = mean(values);
    medianValue = median(values);
    stdValue = std(values);
    ciLow = linear_percentile(values, 0.025);
    ciHigh = linear_percentile(values, 0.975);
end

function percentileValue = linear_percentile(values, probability)
    values = sort(double(values(:)));
    if numel(values) == 1
        percentileValue = values(1);
        return;
    end
    position = 1+(numel(values)-1)*double(probability);
    lowerIndex = floor(position);
    upperIndex = ceil(position);
    fraction = position-lowerIndex;
    percentileValue = values(lowerIndex) + ...
        fraction*(values(upperIndex)-values(lowerIndex));
end

%% Local utility: reproducible seeds
function seed = trial_optimizer_seed(cfg, direction, trialNumber)
    directionOffset = direction_seed_offset(direction);
    seed = valid_seed(double(cfg.optimizationBaseSeed) + ...
        directionOffset + round(double(trialNumber)));
end

function seed = stability_noise_seed( ...
        cfg, direction, trialNumber, levelIndex, replicateIndex)
    directionOffset = direction_seed_offset(direction);
    seed = valid_seed(double(cfg.noiseStabilityBaseSeed) + ...
        directionOffset + 1000*round(double(trialNumber)) + ...
        100000*double(levelIndex) + ...
        10000000*double(replicateIndex));
end

function offset = direction_seed_offset(direction)
    if strcmpi(scalar_string(direction), "Right")
        offset = 0;
    elseif strcmpi(scalar_string(direction), "Left")
        offset = 100000000;
    else
        offset = 200000000;
    end
end

function seed = valid_seed(value)
    maximumSeed = double(intmax('uint32'));
    seed = mod(round(double(value)), maximumSeed);
end

function fraction = safe_noise_fraction(noiseStd, observedOriginal)
    finiteValues = double(observedOriginal(:));
    finiteValues = finiteValues(isfinite(finiteValues));
    signalStd = std(finiteValues);
    if isempty(noiseStd) || ~isfinite(noiseStd) || signalStd <= eps
        fraction = NaN;
    else
        fraction = double(noiseStd)/signalStd;
    end
end

function value = scalar_string(value)
    while iscell(value) && isscalar(value)
        value = value{1};
    end
    value = string(value);
    value = value(:);
    if isempty(value)
        value = "";
    else
        value = value(1);
    end
end

function value = scalar_char(value)
    value = char(scalar_string(value));
end

function value = get_output_value(outputStructure, fieldName, defaultValue)
    value = defaultValue;
    if isstruct(outputStructure) && isfield(outputStructure, fieldName)
        candidate = outputStructure.(fieldName);
        if isnumeric(candidate) && isscalar(candidate) && isfinite(candidate)
            value = double(candidate);
        end
    end
end

%% Empty result structures
function row = empty_summary_row()
    row = struct( ...
        'Dataset', "", ...
        'Channel', "", ...
        'Direction', "", ...
        'TrialNumber', NaN, ...
        'OriginalT1Ms', NaN, ...
        'OriginalNcoup', NaN, ...
        'OriginalR2', NaN, ...
        'FsHz', NaN, ...
        'FitStartMs', NaN, ...
        'FitEndMs', NaN, ...
        'NoiseStdFraction', NaN, ...
        'NoiseStd', NaN, ...
        'NoiseSeed', NaN, ...
        'OptimizerSeed', NaN, ...
        'OptimizationSucceeded', false, ...
        'Converged', false, ...
        'ExitFlag', NaN, ...
        'TotalStarts', 0, ...
        'SuccessfulStarts', 0, ...
        'AllInsideFitWindow', false, ...
        'WindowValid', false, ...
        'ThreeAtomsActive', false, ...
        'ThreeDelayValid', false, ...
        'Delta1Ms', NaN, ...
        'Delta2Ms', NaN, ...
        'Delta3Ms', NaN, ...
        'MinimumSeparationMs', NaN, ...
        'AuxiliaryCoefficient1', NaN, ...
        'AuxiliaryCoefficient2', NaN, ...
        'AuxiliaryCoefficient3', NaN, ...
        'MinimumAuxiliaryShare', NaN, ...
        'Baseline', NaN, ...
        'ProfiledR2', NaN, ...
        'ProfiledRMSE', NaN, ...
        'ProfiledSSE', NaN, ...
        'Status', "", ...
        'Message', "");
end

function model = empty_model_result()
    model = struct( ...
        'numberOfAtoms', 0, ...
        'delays', [], ...
        'delaysMs', [], ...
        'amplitudes', [], ...
        'auxiliaryShares', [], ...
        'baseline', NaN, ...
        'atomMatrix', [], ...
        'fittedSignal', [], ...
        'residual', [], ...
        'SSE', inf, ...
        'R2', NaN, ...
        'RMSE', inf, ...
        'exitFlag', NaN, ...
        'converged', false, ...
        'optimizationSucceeded', false, ...
        'iterations', NaN, ...
        'functionEvaluations', NaN, ...
        'totalStarts', 0, ...
        'successfulStarts', 0, ...
        'totalAllStageStarts', 0, ...
        'successfulAllStageStarts', 0);
end

function row = empty_stability_run_row()
    row = struct( ...
        'Dataset', "", ...
        'Channel', "", ...
        'Direction', "", ...
        'TrialNumber', NaN, ...
        'NoiseLevelFraction', NaN, ...
        'NoiseLevelPercent', NaN, ...
        'Replicate', NaN, ...
        'NoiseStd', NaN, ...
        'NoiseSeed', NaN, ...
        'BaselineThreeDelayValid', false, ...
        'OptimizationSucceeded', false, ...
        'Converged', false, ...
        'ThreeAtomsActive', false, ...
        'ThreeDelayValid', false, ...
        'Delta1Ms', NaN, ...
        'Delta2Ms', NaN, ...
        'Delta3Ms', NaN, ...
        'Delta1AbsoluteChangeMs', NaN, ...
        'Delta2AbsoluteChangeMs', NaN, ...
        'Delta3AbsoluteChangeMs', NaN, ...
        'MaximumAbsoluteDelayChangeMs', NaN, ...
        'DelaysWithinToleranceOfBaseline', false, ...
        'Status', "", ...
        'Message', "");
end

function row = empty_stability_summary_row()
    row = struct( ...
        'Dataset', "", ...
        'Channel', "", ...
        'Direction', "", ...
        'TrialNumber', NaN, ...
        'NoiseLevelFraction', NaN, ...
        'NoiseLevelPercent', NaN, ...
        'ReplicateCount', 0, ...
        'DelayToleranceMs', NaN, ...
        'BaselineThreeDelayValid', false, ...
        'BaselineDelta1Ms', NaN, ...
        'BaselineDelta2Ms', NaN, ...
        'BaselineDelta3Ms', NaN, ...
        'OptimizationSucceededCount', 0, ...
        'OptimizationSucceededProportion', NaN, ...
        'ConvergedCount', 0, ...
        'ConvergedProportion', NaN, ...
        'ThreeAtomsActiveCount', 0, ...
        'ThreeAtomsActiveProportion', NaN, ...
        'ThreeDelayValidCount', 0, ...
        'ThreeDelayValidProportion', NaN, ...
        'DelayStableCount', 0, ...
        'DelayStableProportionAmongValid', NaN, ...
        'DelayStableProportionAllReplicates', NaN, ...
        'MeanMaximumAbsoluteDelayChangeMs', NaN, ...
        'MedianMaximumAbsoluteDelayChangeMs', NaN, ...
        'Delta1MeanMs', NaN, ...
        'Delta1MedianMs', NaN, ...
        'Delta1StdMs', NaN, ...
        'Delta1CI95LowMs', NaN, ...
        'Delta1CI95HighMs', NaN, ...
        'Delta2MeanMs', NaN, ...
        'Delta2MedianMs', NaN, ...
        'Delta2StdMs', NaN, ...
        'Delta2CI95LowMs', NaN, ...
        'Delta2CI95HighMs', NaN, ...
        'Delta3MeanMs', NaN, ...
        'Delta3MedianMs', NaN, ...
        'Delta3StdMs', NaN, ...
        'Delta3CI95LowMs', NaN, ...
        'Delta3CI95HighMs', NaN);
end
