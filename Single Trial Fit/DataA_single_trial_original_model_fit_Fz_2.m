%% Data A: robust single-trial original-model fitting with fixed ranges
% This script applies the final workflow established in Steps 2--4.
%
% Main features
%   1. Fits every retained Right- and Left-cue trial independently.
%   2. Uses the original coupling-release model without tau_d.
%   3. Keeps the cue-locked time axis throughout the analysis.
%   4. Fixes sigma0 = 0.1 and K0 = 2.0 for every trial; optimizes only
%      beta, lambda, T_K, p_K, and t1 in scaled coordinates in [0,1].
%   5. Builds a complete integer N_coup profile with forward and reverse
%      continuation, then re-optimizes only suspicious profile drops.
%   6. Selects the first stable multi-point R-squared plateau, rather than
%      accepting an isolated high point or a negligible upper-edge gain.
%   7. Uses the fixed continuous-parameter ranges requested for the
%      present analysis; no boundary-hit diagnosis or automatic expansion
%      is performed.
%   8. Reports R-squared, RMSE, MAE, correlation, SSE, convergence,
%      and N-profile stability in the Command Window and saved tables.
%   9. Summarizes parameter median, Q1, Q3, and IQR separately for
%      Right and Left cues.
%  10. Produces only three final figures: Right-cue fits, Left-cue fits,
%      and the parameter distributions. Fit panels share axis labels;
%      each panel title contains only the trial number and R-squared.
%
% Required input
%   DataA_Fz_single_trial_preprocessed.mat
%
% The input file must contain the structure DataA produced by
%   DataA_single_trial_preprocessing_Fz.m
%
% Required toolbox
%   Optimization Toolbox (lsqnonlin)

clearvars;
close all;
clc;

%% =========================================================
%  File names
%% =========================================================

inputFile = 'DataA_Fz_single_trial_preprocessed.mat';

outputMatFile = ...
    'DataA_Fz_single_trial_original_model_fit_fixed_ranges.mat';

outputTrialCsvFile = ...
    'DataA_Fz_single_trial_fit_summary_fixed_ranges.csv';

outputParameterCsvFile = ...
    'DataA_Fz_parameter_summary_by_cue_fixed_ranges.csv';

ExportFigures = false;

%% =========================================================
%  Original fixed model settings
%% =========================================================

K1 = 2.2;
gamma = 1;
Kc = 2.0;
FixedSigma0 = 0.1;
FixedK0 = 2.0;

% The old average-waveform fit used t1 = 0.075--0.125 s and N = 4--7.
% In the present cue-locked analysis, t1 is usually about 0.5--1.0 s.
% Because delta_t = t1/N, the old N range is no longer appropriate.
% The pilot uses a step of four. The final search evaluates every integer
% N by bidirectional continuation, so the plotted profile has no artificial
% gaps between coarse candidates.
InitialNCandidates = 4:4:80;

parameterNames = { ...
    'sigma0', ...
    'beta', ...
    'lambda', ...
    'K0', ...
    'T_K', ...
    'p_K', ...
    't1'};

% Full parameter layout is retained for model evaluation and saved tables.
% Equal lower/upper bounds encode fixed sigma0 and K0; these entries are
% excluded from optimizer coordinates. Only beta, lambda, T_K, p_K, and
% trial-specific t1 are fitted (N_coup is searched separately as an integer).
InitialLower6 = [ ...
    FixedSigma0, ...
    0.10, ...
    0.10, ...
    FixedK0, ...
    0.050, ...
    0.50];

InitialUpper6 = [ ...
    FixedSigma0, ...
    5.00, ...
    5.00, ...
    FixedK0, ...
    0.20, ...
    2.50];

% Free-parameter starts follow the original GlobalSearch choices.
% Fixed entries equal their prescribed values in every search start.
Center6 = [ ...
    FixedSigma0, ...
    2.00, ...
    2.00, ...
    FixedK0, ...
    0.100, ...
    1.25];

%% =========================================================
%  Robust profile-search settings
%% =========================================================

% Pilot stage:
% Each N_coup value receives one center start and one random start.
PilotStartsPerN = 2;

% Each integer N is traced once from low to high and once from high to low.
% The two directions provide two independent continuation branches without
% running a large blind multi-start search at every N.
ContinuationStartsPerDirection = 1;

% A point is suspicious when its R-squared is this far below both immediate
% neighbours. Only such points receive extra multi-start rescue.
ProfileDropTolerance = 0.010;
MaxProfileRescuePasses = 2;
ProfileRescueStarts = 5;
TopProfileCandidatesToRefine = 3;
TopProfileExtraStarts = 3;

% When several N values are within this absolute R-squared difference of
% the best profile value, use the smallest one. This prevents tiny numerical
% improvements from forcing N to the largest searched value.
PlateauR2Tolerance = 0.003;
PlateauRunLength = 3;

% If the best R-squared near the upper edge is still more than the plateau
% tolerance above all values at least this far below the upper edge, the
% N range is flagged as inadequate.
NUpperEdgeWidth = 8;

% t1 is cue locked and is allowed over the complete positive fitting
% window. The old GlobalSearch t1 interval 0.075--0.125 s was defined
% after resetting the positive-wave onset to t = 0; it is therefore not
% applicable to the present cue-locked single-trial time axis.
T1HalfWidthSeconds = inf;

MinimumFitSamples = 8;
RandomSeedBase = 1;

%% =========================================================
%  Optimizer settings
%% =========================================================

if exist('lsqnonlin', 'file') ~= 2
    error(['This script requires lsqnonlin from Optimization Toolbox. ', ...
           'The original scripts also require Optimization Toolbox.']);
end

optimizerOptions = optimoptions('lsqnonlin', ...
    'Display', 'off', ...
    'Algorithm', 'trust-region-reflective', ...
    'MaxIterations', 300, ...
    'MaxFunctionEvaluations', 3000, ...
    'FunctionTolerance', 1e-8, ...
    'OptimalityTolerance', 1e-8, ...
    'StepTolerance', 1e-8);

settings = struct;
settings.K1 = K1;
settings.gamma = gamma;
settings.Kc = Kc;
settings.PilotStartsPerN = PilotStartsPerN;
settings.ContinuationStartsPerDirection = ...
    ContinuationStartsPerDirection;
settings.ProfileDropTolerance = ProfileDropTolerance;
settings.MaxProfileRescuePasses = MaxProfileRescuePasses;
settings.ProfileRescueStarts = ProfileRescueStarts;
settings.TopProfileCandidatesToRefine = ...
    TopProfileCandidatesToRefine;
settings.TopProfileExtraStarts = TopProfileExtraStarts;
settings.PlateauR2Tolerance = PlateauR2Tolerance;
settings.PlateauRunLength = PlateauRunLength;
settings.NUpperEdgeWidth = NUpperEdgeWidth;
settings.T1HalfWidthSeconds = T1HalfWidthSeconds;
settings.MinimumFitSamples = MinimumFitSamples;
settings.RandomSeedBase = RandomSeedBase;
settings.OptimizerOptions = optimizerOptions;

%% =========================================================
%  Load and validate the preprocessing result
%% =========================================================

if ~isfile(inputFile)
    [selectedFile, selectedPath] = uigetfile( ...
        '*.mat', ...
        'Select DataA_Fz_single_trial_preprocessed.mat');

    if isequal(selectedFile, 0)
        error('No preprocessing MAT file was selected.');
    end

    inputFile = fullfile(selectedPath, selectedFile);
end

loadedData = load(inputFile);

if ~isfield(loadedData, 'DataA')
    error('The selected MAT file does not contain a structure named DataA.');
end

DataA = loadedData.DataA;
validate_data_structure(DataA);

trials = collect_trials(DataA);
numberOfTrials = numel(trials);

fprintf('\nLoaded %d retained trials: %d Right and %d Left.\n', ...
    numberOfTrials, ...
    sum(strcmp({trials.direction}, 'Right')), ...
    sum(strcmp({trials.direction}, 'Left')));

%% =========================================================
%  Stage 1: inexpensive pilot search
%% =========================================================

fprintf('\n============================================================\n');
fprintf('Stage 1: fast pilot search\n');
fprintf('============================================================\n');
fprintf(['Each trial uses %d starts for each of %d initial ', ...
         'N_coup values.\n'], ...
    PilotStartsPerN, numel(InitialNCandidates));

pilotResults = repmat(empty_trial_result(), numberOfTrials, 1);
stage1Timer = tic;

for i = 1:numberOfTrials

    trialTimer = tic;

    randomSeed = ...
        RandomSeedBase + 10000 + i;

    pilotResults(i) = fit_one_trial( ...
        trials(i), ...
        InitialLower6, ...
        InitialUpper6, ...
        Center6, ...
        InitialNCandidates, ...
        settings, ...
        'pilot', ...
        [], ...
        randomSeed);

    elapsedTrial = toc(trialTimer);
    elapsedStage = toc(stage1Timer);
    remaining = ...
        elapsedStage / i * (numberOfTrials - i);

    print_trial_progress( ...
        'Pilot', i, numberOfTrials, pilotResults(i), ...
        elapsedTrial, remaining);
end

stage1ElapsedSeconds = toc(stage1Timer);

%% =========================================================
%  Fixed final ranges
%% =========================================================

% Continuous ranges are fixed throughout the analysis. The pilot solutions
% are used only as warm starts for the robust full-integer N search.
FinalLower6 = InitialLower6;
FinalUpper6 = InitialUpper6;
FinalNCandidates = InitialNCandidates;

fprintf('\nFixed continuous lower bounds: %s\n', ...
    mat2str(FinalLower6, 6));
fprintf('Fixed continuous upper bounds: %s\n', ...
    mat2str(FinalUpper6, 6));
fprintf('Final N_coup span: %d:%d (every integer in Stage 2).\n', ...
    min(FinalNCandidates), max(FinalNCandidates));

SearchRangeTable = table( ...
    string(parameterNames(1:6)).', ...
    FinalLower6(:), ...
    FinalUpper6(:), ...
    'VariableNames', {'Parameter', 'Lower', 'Upper'});

disp(' ');
disp('Parameter bounds (equal lower/upper values denote fixed parameters):');
disp(SearchRangeTable);

%% =========================================================
%  Stage 2: robust full-integer N profiles
%% =========================================================

fprintf('\n============================================================\n');
fprintf('Stage 2: robust bidirectional N-profile search\n');
fprintf('============================================================\n');
fprintf(['Each trial is traced across every integer N in both ', ...
         'directions. Suspicious profile drops are re-optimized only ', ...
         'when necessary. Continuous parameter ranges remain fixed.\n']);

finalResults = repmat(empty_trial_result(), numberOfTrials, 1);
stage2Timer = tic;

for i = 1:numberOfTrials

    trialTimer = tic;
    randomSeed = RandomSeedBase + 20000 + i;

    finalResults(i) = fit_one_trial( ...
        trials(i), ...
        FinalLower6, ...
        FinalUpper6, ...
        Center6, ...
        FinalNCandidates, ...
        settings, ...
        'final', ...
        pilotResults(i), ...
        randomSeed);

    elapsedTrial = toc(trialTimer);
    elapsedStage = toc(stage2Timer);
    remaining = ...
        elapsedStage / i * (numberOfTrials - i);

    print_trial_progress( ...
        'Final', i, numberOfTrials, finalResults(i), ...
        elapsedTrial, remaining);
end

stage2ElapsedSeconds = toc(stage2Timer);

%% =========================================================
%  Create trial-level and cue-level summaries
%% =========================================================

TrialSummary = create_trial_summary_table(finalResults);

ParameterSummary = create_parameter_summary_table( ...
    finalResults, parameterNames);

disp(' ');
disp('Final trial-level fitting summary:');
disp(TrialSummary);

disp(' ');
disp('Parameter summary by cue direction:');
disp(ParameterSummary);

%% =========================================================
%  Save numerical results
%% =========================================================

FitAnalysis = struct;
FitAnalysis.dataset = get_optional_field(DataA, 'dataset', 'Data A');
FitAnalysis.channel = get_optional_field(DataA, 'channel', 'Fz');
FitAnalysis.description = ...
    ['Independent single-trial fitting with the original coupling-release ', ...
     'model, cue-locked time, fixed sigma0 = 0.1 and K0 = 2.0, ', ...
     'and bounded fitting of beta, lambda, T_K, p_K, and t1. ', ...
     'No tau_d parameter or parameter-boundary diagnosis was included.'];

FitAnalysis.model.sigma0 = FixedSigma0;
FitAnalysis.model.K0 = FixedK0;
FitAnalysis.model.fixed_parameter_names = {'sigma0', 'K0'};
FitAnalysis.model.fitted_parameter_names = {'beta', 'lambda', 'T_K', 'p_K', 't1'};
FitAnalysis.model.K1 = K1;
FitAnalysis.model.gamma = gamma;
FitAnalysis.model.Kc = Kc;
FitAnalysis.model.parameter_names = parameterNames;

FitAnalysis.search.InitialLower6 = InitialLower6;
FitAnalysis.search.InitialUpper6 = InitialUpper6;
FitAnalysis.search.InitialNCandidates = InitialNCandidates;
FitAnalysis.search.FinalLower6 = FinalLower6;
FitAnalysis.search.FinalUpper6 = FinalUpper6;
FitAnalysis.search.FinalNCandidates = FinalNCandidates;
FitAnalysis.search.RangeTable = SearchRangeTable;
FitAnalysis.search.settings = settings;

FitAnalysis.elapsed.PilotSeconds = stage1ElapsedSeconds;
FitAnalysis.elapsed.FinalSeconds = stage2ElapsedSeconds;
FitAnalysis.elapsed.TotalSeconds = ...
    stage1ElapsedSeconds + stage2ElapsedSeconds;

FitAnalysis.pilotResults = pilotResults;
FitAnalysis.finalResults = finalResults;
FitAnalysis.trialSummary = TrialSummary;
FitAnalysis.parameterSummary = ParameterSummary;

save(outputMatFile, 'FitAnalysis', '-v7');
writetable(TrialSummary, outputTrialCsvFile);
writetable(ParameterSummary, outputParameterCsvFile);

fprintf('\n============================================================\n');
fprintf('Analysis completed in %.2f minutes.\n', ...
    FitAnalysis.elapsed.TotalSeconds / 60);
fprintf('Saved complete results: %s\n', outputMatFile);
fprintf('Saved trial summary:    %s\n', outputTrialCsvFile);
fprintf('Saved parameter summary:%s\n', outputParameterCsvFile);

%% =========================================================
%  Plot final results
%% =========================================================

figRight = plot_direction_fits(finalResults, 'Right', 1);
figLeft = plot_direction_fits(finalResults, 'Left', 2);
figParameters = plot_parameter_distributions( ...
    finalResults, parameterNames, 3);

if ExportFigures

    exportgraphics(figRight, ...
        'DataA_right_single_trial_original_model_fits_fixed_ranges.pdf', ...
        'ContentType', 'vector', ...
        'BackgroundColor', 'white');

    exportgraphics(figLeft, ...
        'DataA_left_single_trial_original_model_fits_fixed_ranges.pdf', ...
        'ContentType', 'vector', ...
        'BackgroundColor', 'white');

    exportgraphics(figParameters, ...
        'DataA_parameter_distributions_by_cue_fixed_ranges.pdf', ...
        'ContentType', 'vector', ...
        'BackgroundColor', 'white');
end

%% =========================================================
%  Local function: validate the DataA structure
%% =========================================================

function validate_data_structure(DataA)

    requiredTopFields = {'right', 'left'};

    for i = 1:numel(requiredTopFields)
        if ~isfield(DataA, requiredTopFields{i})
            error('DataA is missing field: %s', requiredTopFields{i});
        end
    end

    requiredGroupFields = { ...
        'trial_number', ...
        'cue_sample', ...
        'peak_amplitude_raw_imf5', ...
        'peak_time_cue_locked', ...
        'fit_start_time_cue_locked', ...
        'fit_end_time_cue_locked', ...
        'fit_time_cue_locked', ...
        'fit_data_normalized'};

    groupNames = {'right', 'left'};

    for g = 1:numel(groupNames)

        group = DataA.(groupNames{g});

        for f = 1:numel(requiredGroupFields)
            if ~isfield(group, requiredGroupFields{f})
                error('%s group is missing field: %s', ...
                    groupNames{g}, requiredGroupFields{f});
            end
        end

        numberOfGroupTrials = numel(group.trial_number);

        if numel(group.fit_time_cue_locked) ~= numberOfGroupTrials || ...
           numel(group.fit_data_normalized) ~= numberOfGroupTrials
            error('%s group contains inconsistent trial counts.', ...
                groupNames{g});
        end
    end
end

%% =========================================================
%  Local function: collect both cue directions
%% =========================================================

function trials = collect_trials(DataA)

    groupNames = {'right', 'left'};
    directionNames = {'Right', 'Left'};

    numberOfTrials = ...
        numel(DataA.right.trial_number) + ...
        numel(DataA.left.trial_number);

    template = struct( ...
        'direction', '', ...
        'trial_number', NaN, ...
        'cue_sample', NaN, ...
        'raw_peak_amplitude', NaN, ...
        'peak_time', NaN, ...
        'fit_start_time', NaN, ...
        'fit_end_time', NaN, ...
        't_data', [], ...
        'y_data', []);

    trials = repmat(template, numberOfTrials, 1);
    trialIndex = 0;

    for g = 1:2

        group = DataA.(groupNames{g});

        for i = 1:numel(group.trial_number)

            trialIndex = trialIndex + 1;

            trials(trialIndex).direction = directionNames{g};
            trials(trialIndex).trial_number = group.trial_number(i);
            trials(trialIndex).cue_sample = group.cue_sample(i);
            trials(trialIndex).raw_peak_amplitude = ...
                group.peak_amplitude_raw_imf5(i);
            trials(trialIndex).peak_time = ...
                group.peak_time_cue_locked(i);
            trials(trialIndex).fit_start_time = ...
                group.fit_start_time_cue_locked(i);
            trials(trialIndex).fit_end_time = ...
                group.fit_end_time_cue_locked(i);
            trials(trialIndex).t_data = ...
                group.fit_time_cue_locked{i}(:);
            trials(trialIndex).y_data = ...
                group.fit_data_normalized{i}(:);
        end
    end
end

%% =========================================================
%  Local function: fit one trial
%% =========================================================

function result = fit_one_trial( ...
    trial, lower6, upper6, center6, NCandidates, settings, ...
    searchMode, pilotResult, randomSeed)

    result = empty_trial_result();

    result.direction = trial.direction;
    result.trial_number = trial.trial_number;
    result.cue_sample = trial.cue_sample;
    result.raw_peak_amplitude = trial.raw_peak_amplitude;
    result.peak_time_cue_locked = trial.peak_time;
    result.fit_start_time_cue_locked = trial.fit_start_time;
    result.fit_end_time_cue_locked = trial.fit_end_time;

    tData = trial.t_data(:);
    yData = trial.y_data(:);

    validMask = isfinite(tData) & isfinite(yData);
    tData = tData(validMask);
    yData = yData(validMask);

    result.t_data = tData;
    result.y_data = yData;
    result.number_of_samples = numel(tData);

    if numel(tData) < settings.MinimumFitSamples
        result.status_message = sprintf( ...
            'Fewer than %d valid samples.', ...
            settings.MinimumFitSamples);
        return
    end

    if any(diff(tData) <= 0)
        result.status_message = ...
            'Cue-locked fitting time is not strictly increasing.';
        return
    end

    if std(yData) <= eps
        result.status_message = 'Fitting data are constant.';
        return
    end

    % t1 remains cue locked and is searched near the observed trial peak.
    t1Lower = max( ...
        trial.fit_start_time, ...
        trial.peak_time - settings.T1HalfWidthSeconds);

    t1Upper = min( ...
        trial.fit_end_time, ...
        trial.peak_time + settings.T1HalfWidthSeconds);

    if ~isfinite(t1Lower) || ~isfinite(t1Upper) || ...
       t1Upper <= t1Lower
        t1Lower = min(tData);
        t1Upper = max(tData);
    end

    lowerBounds = [lower6, t1Lower];
    upperBounds = [upper6, t1Upper];
    center = [center6, trial.peak_time];
    center = clip_to_bounds(center, lowerBounds, upperBounds);

    result.lower_bounds = lowerBounds;
    result.upper_bounds = upperBounds;

    rng(randomSeed, 'twister');

    numberOfN = numel(NCandidates);
    perN = repmat(empty_candidate(), numberOfN, 1);

    if strcmp(searchMode, 'pilot')

        for nIndex = 1:numberOfN

            NCurrent = NCandidates(nIndex);

            startPoints = make_start_points( ...
                settings.PilotStartsPerN, ...
                center, ...
                [], ...
                lowerBounds, ...
                upperBounds);

            perN(nIndex) = search_one_N( ...
                NCurrent, ...
                startPoints, ...
                tData, ...
                yData, ...
                lowerBounds, ...
                upperBounds, ...
                settings);
        end

    elseif strcmp(searchMode, 'final')

        % Build a complete integer profile. A full profile is important:
        % connecting sparse coarse candidates can conceal failed local
        % optimizations and can make N selection depend on grid spacing.
        profileN = min(NCandidates):max(NCandidates);
        numberOfN = numel(profileN);
        perN = repmat(empty_candidate(), numberOfN, 1);

        % Phase A: low-to-high continuation.
        for nIndex = 1:numberOfN

            NCurrent = profileN(nIndex);
            neighbourStart = [];

            if nIndex > 1 && perN(nIndex - 1).success
                neighbourStart = perN(nIndex - 1).parameters;
            end

            pilotStart = ...
                get_pilot_warm_start(pilotResult, NCurrent);

            if isempty(pilotStart) && ...
               ~isempty(pilotResult) && ...
               isfield(pilotResult, 'per_N')
                pilotStart = nearest_candidate_start( ...
                    pilotResult.per_N, NCurrent);
            end

            % At exact pilot N values, re-optimize the independently found
            % pilot solution. Between pilot points, continue from N-1.
            if ~isempty(get_pilot_warm_start(pilotResult, NCurrent))
                primaryStart = pilotStart;
                secondaryStart = neighbourStart;
            else
                primaryStart = neighbourStart;
                secondaryStart = pilotStart;
            end

            startPoints = make_continuation_start_points( ...
                settings.ContinuationStartsPerDirection, ...
                primaryStart, ...
                secondaryStart, ...
                center, ...
                lowerBounds, ...
                upperBounds);

            perN(nIndex) = search_one_N( ...
                NCurrent, ...
                startPoints, ...
                tData, ...
                yData, ...
                lowerBounds, ...
                upperBounds, ...
                settings);
        end

        % Phase B: high-to-low continuation. The candidate already found
        % at the same N is retained as one seed, while the next-higher N
        % supplies an independent continuation branch.
        for nIndex = numberOfN:-1:1

            NCurrent = profileN(nIndex);
            primaryStart = [];

            if nIndex < numberOfN && perN(nIndex + 1).success
                primaryStart = perN(nIndex + 1).parameters;
            end

            secondaryStart = [];
            if perN(nIndex).success
                secondaryStart = perN(nIndex).parameters;
            end

            startPoints = make_continuation_start_points( ...
                settings.ContinuationStartsPerDirection, ...
                primaryStart, ...
                secondaryStart, ...
                center, ...
                lowerBounds, ...
                upperBounds);

            reverseCandidate = search_one_N( ...
                NCurrent, ...
                startPoints, ...
                tData, ...
                yData, ...
                lowerBounds, ...
                upperBounds, ...
                settings);

            perN(nIndex) = choose_better_candidate( ...
                perN(nIndex), reverseCandidate);
        end

        % Phase C: identify abrupt downward jumps. These points are not
        % smoothed or replaced numerically; they are genuinely re-fitted
        % using solutions from both neighbours and the current global best.
        [anomalyIndices, ~] = find_profile_anomalies( ...
            perN, settings.ProfileDropTolerance);
        result.profile_jump_count_before = numel(anomalyIndices);

        for rescuePass = 1:settings.MaxProfileRescuePasses

            if isempty(anomalyIndices)
                break
            end

            successfulR2 = [perN.R2];
            successfulR2(~[perN.success]) = -inf;
            [~, globalBestIndex] = max(successfulR2);

            for a = 1:numel(anomalyIndices)

                nIndex = anomalyIndices(a);
                seedRows = zeros(0, numel(lowerBounds));

                if perN(nIndex).success
                    seedRows(end + 1, :) = ...
                        perN(nIndex).parameters; %#ok<AGROW>
                end

                if nIndex > 1 && perN(nIndex - 1).success
                    seedRows(end + 1, :) = ...
                        perN(nIndex - 1).parameters; %#ok<AGROW>
                end

                if nIndex < numberOfN && perN(nIndex + 1).success
                    seedRows(end + 1, :) = ...
                        perN(nIndex + 1).parameters; %#ok<AGROW>
                end

                if isfinite(globalBestIndex) && ...
                   perN(globalBestIndex).success
                    seedRows(end + 1, :) = ...
                        perN(globalBestIndex).parameters; %#ok<AGROW>
                end

                pilotStart = get_pilot_warm_start( ...
                    pilotResult, profileN(nIndex));
                if ~isempty(pilotStart)
                    seedRows(end + 1, :) = pilotStart; %#ok<AGROW>
                end

                startPoints = make_multi_seed_start_points( ...
                    settings.ProfileRescueStarts, ...
                    seedRows, ...
                    center, ...
                    lowerBounds, ...
                    upperBounds);

                rescuedCandidate = search_one_N( ...
                    profileN(nIndex), ...
                    startPoints, ...
                    tData, ...
                    yData, ...
                    lowerBounds, ...
                    upperBounds, ...
                    settings);

                previousCandidate = perN(nIndex);
                perN(nIndex) = choose_better_candidate( ...
                    previousCandidate, rescuedCandidate);

                if perN(nIndex).R2 > previousCandidate.R2 + 1e-10
                    perN(nIndex).rescue_count = ...
                        previousCandidate.rescue_count + 1;
                end
            end

            [anomalyIndices, ~] = find_profile_anomalies( ...
                perN, settings.ProfileDropTolerance);
        end

        % Phase D: give the best few N values independent local and random
        % starts. This protects the plateau height from a shared continuation
        % branch that is smooth but suboptimal.
        r2ByN = [perN.R2];
        r2ByN(~[perN.success]) = -inf;
        [~, ranking] = sort(r2ByN, 'descend');
        numberToRefine = min( ...
            settings.TopProfileCandidatesToRefine, ...
            sum(isfinite(r2ByN)));

        for rankIndex = 1:numberToRefine

            nIndex = ranking(rankIndex);
            startPoints = make_start_points( ...
                settings.TopProfileExtraStarts, ...
                center, ...
                perN(nIndex).parameters, ...
                lowerBounds, ...
                upperBounds);

            refinedCandidate = search_one_N( ...
                profileN(nIndex), ...
                startPoints, ...
                tData, ...
                yData, ...
                lowerBounds, ...
                upperBounds, ...
                settings);

            perN(nIndex) = choose_better_candidate( ...
                perN(nIndex), refinedCandidate);
        end

        [remainingAnomalies, ~] = find_profile_anomalies( ...
            perN, settings.ProfileDropTolerance);

        % A top-candidate refinement can occasionally raise one point enough
        % to expose a neighbouring notch. Give those newly exposed points one
        % last deterministic neighbour-based rescue.
        if ~isempty(remainingAnomalies)

            r2ByN = [perN.R2];
            r2ByN(~[perN.success]) = -inf;
            [~, globalBestIndex] = max(r2ByN);

            for a = 1:numel(remainingAnomalies)

                nIndex = remainingAnomalies(a);
                seedRows = zeros(0, numel(lowerBounds));

                if perN(nIndex).success
                    seedRows(end + 1, :) = ...
                        perN(nIndex).parameters; %#ok<AGROW>
                end

                if nIndex > 1 && perN(nIndex - 1).success
                    seedRows(end + 1, :) = ...
                        perN(nIndex - 1).parameters; %#ok<AGROW>
                end

                if nIndex < numberOfN && perN(nIndex + 1).success
                    seedRows(end + 1, :) = ...
                        perN(nIndex + 1).parameters; %#ok<AGROW>
                end

                if perN(globalBestIndex).success
                    seedRows(end + 1, :) = ...
                        perN(globalBestIndex).parameters; %#ok<AGROW>
                end

                startPoints = make_multi_seed_start_points( ...
                    settings.ProfileRescueStarts, ...
                    seedRows, ...
                    center, ...
                    lowerBounds, ...
                    upperBounds);

                rescuedCandidate = search_one_N( ...
                    profileN(nIndex), ...
                    startPoints, ...
                    tData, ...
                    yData, ...
                    lowerBounds, ...
                    upperBounds, ...
                    settings);

                previousCandidate = perN(nIndex);
                perN(nIndex) = choose_better_candidate( ...
                    previousCandidate, rescuedCandidate);

                if perN(nIndex).R2 > previousCandidate.R2 + 1e-10
                    perN(nIndex).rescue_count = ...
                        previousCandidate.rescue_count + 1;
                end
            end
        end

        [remainingAnomalies, ~] = find_profile_anomalies( ...
            perN, settings.ProfileDropTolerance);
        result.profile_jump_count_after = ...
            numel(remainingAnomalies);
        result.N_profile_stable = isempty(remainingAnomalies);

    else
        error('Unknown search mode: %s', searchMode);
    end

    result.per_N = perN;

    [bestCandidate, plateauReferenceR2] = select_plateau_candidate( ...
        perN, ...
        settings.PlateauR2Tolerance, ...
        settings.PlateauRunLength);

    if ~bestCandidate.success
        result.status_message = ...
            'No valid optimizer solution was found.';
        return
    end

    result.success = true;
    result.status_message = 'Success';
    result.parameters = bestCandidate.parameters;
    result.N_coup = bestCandidate.N_coup;
    result.y_model = bestCandidate.y_model;
    result.residual = yData - bestCandidate.y_model;
    result.metrics = bestCandidate.metrics;
    result.exitflag = bestCandidate.exitflag;
    result.converged = bestCandidate.converged;
    result.iterations = bestCandidate.iterations;
    result.function_evaluations = ...
        bestCandidate.function_evaluations;
    result.first_order_optimality = ...
        bestCandidate.first_order_optimality;

    result.total_starts = sum([perN.total_starts]);
    result.valid_starts = sum([perN.valid_starts]);
    result.converged_starts = sum([perN.converged_starts]);

    result.N_at_lower_candidate = ...
        result.N_coup == min([perN.N_coup]);

    result.N_at_upper_candidate = ...
        result.N_coup == max([perN.N_coup]);

    successfulMask = [perN.success];
    successfulCandidates = perN(successfulMask);

    if ~isempty(successfulCandidates)
        successfulR2 = [successfulCandidates.R2];
        [result.max_R2_across_N, maxIndex] = max(successfulR2);
        result.N_of_max_R2 = ...
            successfulCandidates(maxIndex).N_coup;
        result.R2_loss_from_max = ...
            result.max_R2_across_N - result.metrics.R2;
    end

    result.N_range_adequate = assess_N_range( ...
        perN, ...
        settings.PlateauR2Tolerance, ...
        settings.NUpperEdgeWidth);
    result.plateau_R2_tolerance = ...
        settings.PlateauR2Tolerance;
    result.plateau_reference_R2 = ...
        plateauReferenceR2;
    result.profile_drop_tolerance = ...
        settings.ProfileDropTolerance;
end

%% =========================================================
%  Local function: search one fixed N_coup value
%% =========================================================

function best = search_one_N( ...
    NCurrent, startPoints, tData, yData, ...
    lowerBounds, upperBounds, settings)

    if isvector(startPoints)
        startPoints = startPoints(:).';
    end

    best = empty_candidate();
    best.N_coup = NCurrent;
    best.total_starts = size(startPoints, 1);

    validStarts = 0;
    convergedStarts = 0;

    for s = 1:size(startPoints, 1)

        x0 = clip_to_bounds( ...
            startPoints(s, :), lowerBounds, upperBounds);

        % Scale only free parameters; zero-width bounds encode constants.
        % This avoids division by zero and keeps sigma0 and K0 out of
        % all pilot, continuation, rescue, and refinement optimizations.
        intervalWidth = upperBounds - lowerBounds;
        freeMask = intervalWidth > 0;
        z0 = (x0(freeMask) - lowerBounds(freeMask)) ./ ...
            intervalWidth(freeMask);
        z0 = min(max(z0, 0), 1);

        residualFunction = @(z) scaled_model_residual( ...
            z, lowerBounds, upperBounds, ...
            NCurrent, tData, yData, ...
            settings.K1, settings.gamma);

        try
            [zFit, ~, ~, exitflag, output] = lsqnonlin( ...
                residualFunction, ...
                z0, ...
                zeros(size(z0)), ...
                ones(size(z0)), ...
                settings.OptimizerOptions);
        catch
            continue
        end

        % Reconstruct the full seven-entry vector for existing outputs.
        xFit = lowerBounds;
        xFit(freeMask) = lowerBounds(freeMask) + ...
            zFit(:).' .* intervalWidth(freeMask);
        xFit = xFit(:).';

        [yModel, metrics, isValid] = evaluate_candidate( ...
            xFit, NCurrent, tData, yData, ...
            settings.K1, settings.gamma);

        if ~isValid
            continue
        end

        validStarts = validStarts + 1;

        if exitflag > 0
            convergedStarts = convergedStarts + 1;
        end

        candidate = empty_candidate();
        candidate.success = true;
        candidate.N_coup = NCurrent;
        candidate.parameters = xFit;
        candidate.y_model = yModel;
        candidate.metrics = metrics;
        candidate.R2 = metrics.R2;
        candidate.RMSE = metrics.RMSE;
        candidate.exitflag = exitflag;
        candidate.converged = exitflag > 0;
        candidate.iterations = get_output_value( ...
            output, 'iterations', NaN);
        candidate.function_evaluations = get_output_value( ...
            output, 'funcCount', NaN);
        candidate.first_order_optimality = get_output_value( ...
            output, 'firstorderopt', NaN);

        best = choose_better_candidate(best, candidate);
    end

    best.total_starts = size(startPoints, 1);
    best.valid_starts = validStarts;
    best.converged_starts = convergedStarts;
end

%% =========================================================
%  Local function: scaled residual vector for lsqnonlin
%% =========================================================

function residual = scaled_model_residual( ...
    z, lowerBounds, upperBounds, ...
    NCurrent, tData, yData, K1, gamma)

    % Fixed entries are copied exactly; only free entries depend on z.
    z = z(:).';
    intervalWidth = upperBounds - lowerBounds;
    freeMask = intervalWidth > 0;
    x = lowerBounds;
    x(freeMask) = lowerBounds(freeMask) + ...
        z .* intervalWidth(freeMask);

    [yModel, ~] = generate_model_curve( ...
        x, NCurrent, tData, K1, gamma);

    if any(~isfinite(yModel)) || ...
       any(imag(yModel) ~= 0) || ...
       max(real(yModel)) > 1.5 || ...
       min(real(yModel)) < -0.2

        residual = 1e3 * ones(size(yData));
        return
    end

    residual = yData(:) - real(yModel(:));
end

%% =========================================================
%  Local function: evaluate one fitted candidate
%% =========================================================

function [yModel, metrics, isValid] = evaluate_candidate( ...
    x, NCurrent, tData, yData, K1, gamma)

    [yModel, ~] = generate_model_curve( ...
        x, NCurrent, tData, K1, gamma);

    isValid = ...
        all(isfinite(x)) && ...
        all(isfinite(yModel)) && ...
        all(imag(yModel) == 0) && ...
        max(real(yModel)) <= 1.5 && ...
        min(real(yModel)) >= -0.2;

    if ~isValid
        yModel = NaN(size(yData));
        metrics = empty_metrics();
        return
    end

    yModel = real(yModel(:));
    errorVector = yData(:) - yModel;

    SSE = sum(errorVector .^ 2);
    SST = sum((yData(:) - mean(yData(:))) .^ 2);

    if SST > 0
        R2 = 1 - SSE / SST;
    else
        R2 = NaN;
    end

    RMSE = sqrt(mean(errorVector .^ 2));
    MAE = mean(abs(errorVector));

    if std(yData) > 0 && std(yModel) > 0
        correlationMatrix = corrcoef(yData(:), yModel);
        Correlation = correlationMatrix(1, 2);
    else
        Correlation = NaN;
    end

    metrics = struct;
    metrics.R2 = R2;
    metrics.RMSE = RMSE;
    metrics.MAE = MAE;
    metrics.Correlation = Correlation;
    metrics.SSE = SSE;

    isValid = ...
        isfinite(R2) && ...
        isfinite(RMSE) && ...
        isfinite(SSE);
end

%% =========================================================
%  Local function: original coupling-release model
%% =========================================================

function [yEvaluation, detail] = generate_model_curve( ...
    x, N_coup, tEvaluation, K1, gamma)

    sigma0 = x(1);
    beta = x(2);
    lambda = x(3);
    K0 = x(4);
    T_K = x(5);
    p_K = x(6);
    t1 = x(7);

    SInfinity = 2/pi;

    tEvaluation = tEvaluation(:);
    tEnd = max(tEvaluation);

    if sigma0 <= 0 || ...
       beta <= 0 || ...
       lambda <= 0 || ...
       K0 <= 0 || ...
       K0 >= K1 || ...
       T_K <= 0 || ...
       p_K <= 0 || ...
       t1 <= 0 || ...
       N_coup < 1

        yEvaluation = NaN(size(tEvaluation));
        detail = struct;
        return
    end

    deltaT = t1 / N_coup;
    tau = deltaT;

    logisticSigma = @(t) SInfinity ./ ...
        (1 + (SInfinity / sigma0 - 1) .* ...
        exp( ...
        -(beta * (8 * K1 / (3 * pi^2 * gamma))) .* ...
        ((1 + t ./ (lambda * tau)) .^ (3/2) - 1) + ...
        t ./ (lambda * tau)));

    tCoupling = (0:N_coup) * deltaT;
    sigmaCoupling = logisticSigma(tCoupling);
    sigmaK = sigmaCoupling(end);

    releaseSteps = max( ...
        3, ...
        ceil((tEnd - t1) / deltaT) + 5);

    tRelease = t1 + (0:releaseSteps) * deltaT;

    KRelease = ...
        K0 + (K1 - K0) .* ...
        exp(-((max(tRelease - t1, 0)) ./ T_K) .^ p_K);

    sigmaRelease = zeros(1, releaseSteps + 1);
    sigmaRelease(1) = sigmaK;

    for j = 1:releaseSteps

        combinedRatio = ...
            (KRelease(j + 1) * sigmaRelease(j)) / ...
            (K1 * sigmaK);

        combinedRatio = max(0, min(combinedRatio, 1));

        sigmaRelease(j + 1) = ...
            (2 * sigmaK / pi) * asin(combinedRatio);
    end

    tAll = [tCoupling, tRelease(2:end)];
    sigmaAll = [sigmaCoupling, sigmaRelease(2:end)];

    maximumSigma = max(sigmaAll);

    if ~isfinite(maximumSigma) || maximumSigma <= 0
        yEvaluation = NaN(size(tEvaluation));
        detail = struct;
        return
    end

    sigmaAll = sigmaAll / maximumSigma;

    yEvaluation = interp1( ...
        tAll, sigmaAll, tEvaluation, 'pchip', 'extrap');

    yEvaluation = yEvaluation(:);
    yEvaluation(yEvaluation < 0) = 0;

    detail = struct;
    detail.t_all = tAll;
    detail.sigma_all = sigmaAll;
    detail.t_coupling = tCoupling;
    detail.sigma_coupling = sigmaCoupling;
    detail.t_release = tRelease;
    detail.sigma_release = sigmaRelease;
    detail.K_release = KRelease;
    detail.delta_t = deltaT;
    detail.t1 = t1;
    detail.sigma_k = sigmaK;
end

%% =========================================================
%  Local function: build search starts
%% =========================================================

function startPoints = make_start_points( ...
    numberOfStarts, center, warmStart, lowerBounds, upperBounds)

    numberOfParameters = numel(lowerBounds);
    startPoints = zeros(numberOfStarts, numberOfParameters);

    if numberOfStarts < 1
        return
    end

    rowIndex = 1;

    if ~isempty(warmStart) && ...
       numel(warmStart) == numberOfParameters && ...
       all(isfinite(warmStart))
        startPoints(rowIndex, :) = ...
            clip_to_bounds(warmStart, lowerBounds, upperBounds);
    else
        startPoints(rowIndex, :) = ...
            clip_to_bounds(center, lowerBounds, upperBounds);
    end

    while rowIndex < numberOfStarts

        rowIndex = rowIndex + 1;

        if ~isempty(warmStart) && mod(rowIndex, 2) == 0

            localStart = warmStart + ...
                0.20 .* (upperBounds - lowerBounds) .* ...
                randn(size(lowerBounds));

            startPoints(rowIndex, :) = ...
                clip_to_bounds( ...
                    localStart, lowerBounds, upperBounds);

        else

            startPoints(rowIndex, :) = ...
                lowerBounds + ...
                rand(size(lowerBounds)) .* ...
                (upperBounds - lowerBounds);
        end
    end
end

%% =========================================================
%  Local function: deterministic continuation starts
%% =========================================================

function startPoints = make_continuation_start_points( ...
    numberOfStarts, primaryStart, secondaryStart, center, ...
    lowerBounds, upperBounds)

    seedRows = zeros(0, numel(lowerBounds));

    if is_valid_start(primaryStart, numel(lowerBounds))
        seedRows(end + 1, :) = primaryStart(:).'; %#ok<AGROW>
    end

    if is_valid_start(secondaryStart, numel(lowerBounds))
        seedRows(end + 1, :) = secondaryStart(:).'; %#ok<AGROW>
    end

    startPoints = make_multi_seed_start_points( ...
        numberOfStarts, ...
        seedRows, ...
        center, ...
        lowerBounds, ...
        upperBounds);
end

%% =========================================================
%  Local function: combine several warm starts without duplicates
%% =========================================================

function startPoints = make_multi_seed_start_points( ...
    numberOfStarts, seedRows, center, lowerBounds, upperBounds)

    numberOfParameters = numel(lowerBounds);

    if isempty(seedRows)
        seedRows = zeros(0, numberOfParameters);
    end

    if size(seedRows, 2) ~= numberOfParameters
        seedRows = zeros(0, numberOfParameters);
    end

    validRows = ...
        all(isfinite(seedRows), 2);
    seedRows = seedRows(validRows, :);

    for i = 1:size(seedRows, 1)
        seedRows(i, :) = clip_to_bounds( ...
            seedRows(i, :), lowerBounds, upperBounds);
    end

    if ~isempty(seedRows)
        % Compare starts only in free coordinates (fixed widths are zero).
        intervalWidth = upperBounds - lowerBounds;
        freeMask = intervalWidth > 0;
        scaledSeeds = ...
            (seedRows(:, freeMask) - lowerBounds(freeMask)) ./ ...
            intervalWidth(freeMask);
        [~, uniqueIndices] = unique( ...
            round(scaledSeeds * 1e10) / 1e10, ...
            'rows', ...
            'stable');
        seedRows = seedRows(sort(uniqueIndices), :);
    end

    center = clip_to_bounds(center, lowerBounds, upperBounds);

    if isempty(seedRows)
        seedRows = center;
    elseif size(seedRows, 1) < numberOfStarts
        seedRows(end + 1, :) = center;
    end

    startPoints = zeros(numberOfStarts, numberOfParameters);
    numberFromSeeds = min(numberOfStarts, size(seedRows, 1));
    startPoints(1:numberFromSeeds, :) = ...
        seedRows(1:numberFromSeeds, :);

    rowIndex = numberFromSeeds;
    while rowIndex < numberOfStarts

        rowIndex = rowIndex + 1;
        referenceRow = ...
            startPoints(mod(rowIndex - 2, numberFromSeeds) + 1, :);

        if mod(rowIndex, 2) == 0
            newStart = referenceRow + ...
                0.12 .* (upperBounds - lowerBounds) .* ...
                randn(1, numberOfParameters);
        else
            newStart = lowerBounds + ...
                rand(1, numberOfParameters) .* ...
                (upperBounds - lowerBounds);
        end

        startPoints(rowIndex, :) = ...
            clip_to_bounds(newStart, lowerBounds, upperBounds);
    end
end

function valid = is_valid_start(startPoint, numberOfParameters)

    valid = ...
        ~isempty(startPoint) && ...
        numel(startPoint) == numberOfParameters && ...
        all(isfinite(startPoint(:)));
end

%% =========================================================
%  Local function: retrieve a pilot warm start for one N
%% =========================================================

function warmStart = get_pilot_warm_start(pilotResult, NCurrent)

    warmStart = [];

    if isempty(pilotResult) || ...
       ~isfield(pilotResult, 'per_N') || ...
       isempty(pilotResult.per_N)
        return
    end

    for i = 1:numel(pilotResult.per_N)

        candidate = pilotResult.per_N(i);

        if candidate.success && candidate.N_coup == NCurrent
            warmStart = candidate.parameters;
            return
        end
    end
end

%% =========================================================
%  Local function: nearest successful candidate warm start
%% =========================================================

function warmStart = nearest_candidate_start(candidates, NCurrent)

    warmStart = [];

    if isempty(candidates)
        return
    end

    successfulMask = [candidates.success];
    successfulCandidates = candidates(successfulMask);

    if isempty(successfulCandidates)
        return
    end

    candidateN = [successfulCandidates.N_coup];
    [~, nearestIndex] = min(abs(candidateN - NCurrent));
    warmStart = successfulCandidates(nearestIndex).parameters;
end

%% =========================================================
%  Local function: select the smallest N on the R-squared plateau
%% =========================================================

function [selected, plateauReferenceR2] = select_plateau_candidate( ...
    candidates, tolerance, runLength)

    selected = empty_candidate();
    plateauReferenceR2 = NaN;

    if isempty(candidates)
        return
    end

    successfulMask = [candidates.success];
    successfulCandidates = candidates(successfulMask);

    if isempty(successfulCandidates)
        return
    end

    r2Values = [successfulCandidates.R2];
    maximumR2 = max(r2Values);

    [NValues, order] = sort([successfulCandidates.N_coup]);
    successfulCandidates = successfulCandidates(order);
    r2Values = [successfulCandidates.R2];

    selectedIndex = [];
    windowStarts = [];
    windowScores = [];

    % Use the median of each consecutive window as the plateau reference.
    % An isolated optimistic point therefore cannot raise the threshold for
    % the whole profile.
    if numel(NValues) >= runLength
        for i = 1:(numel(NValues) - runLength + 1)
            runIndices = i:(i + runLength - 1);
            if all(diff(NValues(runIndices)) == 1)
                windowStarts(end + 1) = i; %#ok<AGROW>
                windowScores(end + 1) = ... %#ok<AGROW>
                    median(r2Values(runIndices));
            end
        end
    end

    if ~isempty(windowScores)
        plateauReferenceR2 = max(windowScores);
        plateauThreshold = ...
            plateauReferenceR2 - tolerance;

        for w = 1:numel(windowStarts)
            runIndices = ...
                windowStarts(w):(windowStarts(w) + runLength - 1);

            if all(r2Values(runIndices) >= plateauThreshold)
                selectedIndex = runIndices(1);
                break
            end
        end
    end

    if isempty(selectedIndex)
        [~, selectedIndex] = max(r2Values);
        plateauReferenceR2 = maximumR2;
    end

    selected = successfulCandidates(selectedIndex);
end

%% =========================================================
%  Local function: find failed or abruptly low N-profile points
%% =========================================================

function [indices, jumpMagnitude] = ...
    find_profile_anomalies(candidates, dropTolerance)

    numberOfCandidates = numel(candidates);
    flagged = false(1, numberOfCandidates);
    jumpMagnitude = zeros(1, numberOfCandidates);

    if numberOfCandidates == 0
        indices = [];
        return
    end

    successMask = [candidates.success];
    r2Values = [candidates.R2];

    flagged(~successMask) = true;

    % Detect a local downward notch relative to both neighbours.
    for i = 2:(numberOfCandidates - 1)
        if successMask(i - 1) && ...
           successMask(i) && ...
           successMask(i + 1)

            localReference = min( ...
                r2Values(i - 1), r2Values(i + 1));
            localDrop = localReference - r2Values(i);

            if localDrop > dropTolerance
                flagged(i) = true;
                jumpMagnitude(i) = ...
                    max(jumpMagnitude(i), localDrop);
            end
        end
    end

    indices = find(flagged);
end

%% =========================================================
%  Local function: check whether the searched N range is adequate
%% =========================================================

function adequate = assess_N_range( ...
    candidates, tolerance, upperEdgeWidth)

    adequate = false;

    if isempty(candidates)
        return
    end

    successfulMask = [candidates.success];
    candidates = candidates(successfulMask);

    if numel(candidates) < 2
        return
    end

    [NValues, order] = sort([candidates.N_coup]);
    candidates = candidates(order);
    r2Values = [candidates.R2];

    maximumN = max(NValues);
    interiorMask = NValues <= maximumN - upperEdgeWidth;

    if ~any(interiorMask)
        return
    end

    bestOverall = max(r2Values);
    bestInterior = max(r2Values(interiorMask));
    [~, bestIndex] = max(r2Values);
    bestN = NValues(bestIndex);

    % The range is adequate when either the best N is clearly away from
    % the upper edge or the improvement near the upper edge is already
    % within the predefined R-squared plateau tolerance.
    adequate = ...
        bestN <= maximumN - upperEdgeWidth || ...
        bestOverall - bestInterior <= tolerance;
end

%% =========================================================
%  Local function: choose the better of two candidates
%% =========================================================

function better = choose_better_candidate(first, second)

    if ~first.success
        better = second;
        return
    end

    if ~second.success
        better = first;
        return
    end

    if second.R2 > first.R2 + 1e-12 || ...
       (abs(second.R2 - first.R2) <= 1e-12 && ...
        second.RMSE < first.RMSE)
        better = second;
    else
        better = first;
    end

    better.total_starts = ...
        first.total_starts + second.total_starts;
    better.valid_starts = ...
        first.valid_starts + second.valid_starts;
    better.converged_starts = ...
        first.converged_starts + second.converged_starts;
    better.rescue_count = ...
        max(first.rescue_count, second.rescue_count);
end

%% =========================================================
%  Local function: clip one point to finite bounds
%% =========================================================

function x = clip_to_bounds(x, lowerBounds, upperBounds)

    x = x(:).';
    x = min(max(x, lowerBounds), upperBounds);
end

%% =========================================================
%  Local function: trial-level summary table
%% =========================================================

function trialTable = create_trial_summary_table(results)

    numberOfTrials = numel(results);

    Direction = strings(numberOfTrials, 1);
    TrialNumber = NaN(numberOfTrials, 1);
    FitSampleCount = NaN(numberOfTrials, 1);
    RawPeakAmplitudeIMF5 = NaN(numberOfTrials, 1);
    PeakTimeCueLocked_s = NaN(numberOfTrials, 1);
    FitStartCueLocked_s = NaN(numberOfTrials, 1);
    FitEndCueLocked_s = NaN(numberOfTrials, 1);

    sigma0 = NaN(numberOfTrials, 1);
    beta = NaN(numberOfTrials, 1);
    lambda = NaN(numberOfTrials, 1);
    K0 = NaN(numberOfTrials, 1);
    T_K = NaN(numberOfTrials, 1);
    p_K = NaN(numberOfTrials, 1);
    t1_s = NaN(numberOfTrials, 1);
    N_coup = NaN(numberOfTrials, 1);
    N_of_max_R2 = NaN(numberOfTrials, 1);
    Max_R2_across_N = NaN(numberOfTrials, 1);
    PlateauReference_R2 = NaN(numberOfTrials, 1);
    R2_loss_from_max = NaN(numberOfTrials, 1);
    NRangeAdequate = false(numberOfTrials, 1);
    NProfileStable = false(numberOfTrials, 1);
    ProfileJumpCountBefore = NaN(numberOfTrials, 1);
    ProfileJumpCountAfter = NaN(numberOfTrials, 1);
    delta_t_s = NaN(numberOfTrials, 1);

    R2 = NaN(numberOfTrials, 1);
    RMSE = NaN(numberOfTrials, 1);
    MAE = NaN(numberOfTrials, 1);
    Correlation = NaN(numberOfTrials, 1);
    SSE = NaN(numberOfTrials, 1);

    Success = false(numberOfTrials, 1);
    Converged = false(numberOfTrials, 1);
    ExitFlag = NaN(numberOfTrials, 1);
    Iterations = NaN(numberOfTrials, 1);
    FunctionEvaluations = NaN(numberOfTrials, 1);
    TotalStarts = NaN(numberOfTrials, 1);
    ValidStarts = NaN(numberOfTrials, 1);
    ConvergedStarts = NaN(numberOfTrials, 1);
    Status = strings(numberOfTrials, 1);

    for i = 1:numberOfTrials

        oneResult = results(i);

        Direction(i) = string(oneResult.direction);
        TrialNumber(i) = oneResult.trial_number;
        FitSampleCount(i) = oneResult.number_of_samples;
        RawPeakAmplitudeIMF5(i) = oneResult.raw_peak_amplitude;
        PeakTimeCueLocked_s(i) = ...
            oneResult.peak_time_cue_locked;
        FitStartCueLocked_s(i) = ...
            oneResult.fit_start_time_cue_locked;
        FitEndCueLocked_s(i) = ...
            oneResult.fit_end_time_cue_locked;
        Success(i) = oneResult.success;
        Status(i) = string(oneResult.status_message);

        if ~oneResult.success
            continue
        end

        sigma0(i) = oneResult.parameters(1);
        beta(i) = oneResult.parameters(2);
        lambda(i) = oneResult.parameters(3);
        K0(i) = oneResult.parameters(4);
        T_K(i) = oneResult.parameters(5);
        p_K(i) = oneResult.parameters(6);
        t1_s(i) = oneResult.parameters(7);
        N_coup(i) = oneResult.N_coup;
        N_of_max_R2(i) = oneResult.N_of_max_R2;
        Max_R2_across_N(i) = oneResult.max_R2_across_N;
        PlateauReference_R2(i) = ...
            oneResult.plateau_reference_R2;
        R2_loss_from_max(i) = oneResult.R2_loss_from_max;
        NRangeAdequate(i) = oneResult.N_range_adequate;
        NProfileStable(i) = oneResult.N_profile_stable;
        ProfileJumpCountBefore(i) = ...
            oneResult.profile_jump_count_before;
        ProfileJumpCountAfter(i) = ...
            oneResult.profile_jump_count_after;
        delta_t_s(i) = t1_s(i) / N_coup(i);

        R2(i) = oneResult.metrics.R2;
        RMSE(i) = oneResult.metrics.RMSE;
        MAE(i) = oneResult.metrics.MAE;
        Correlation(i) = oneResult.metrics.Correlation;
        SSE(i) = oneResult.metrics.SSE;

        Converged(i) = oneResult.converged;
        ExitFlag(i) = oneResult.exitflag;
        Iterations(i) = oneResult.iterations;
        FunctionEvaluations(i) = ...
            oneResult.function_evaluations;
        TotalStarts(i) = oneResult.total_starts;
        ValidStarts(i) = oneResult.valid_starts;
        ConvergedStarts(i) = oneResult.converged_starts;
    end

    trialTable = table( ...
        Direction, ...
        TrialNumber, ...
        FitSampleCount, ...
        RawPeakAmplitudeIMF5, ...
        PeakTimeCueLocked_s, ...
        FitStartCueLocked_s, ...
        FitEndCueLocked_s, ...
        sigma0, ...
        beta, ...
        lambda, ...
        K0, ...
        T_K, ...
        p_K, ...
        t1_s, ...
        N_coup, ...
        N_of_max_R2, ...
        Max_R2_across_N, ...
        PlateauReference_R2, ...
        R2_loss_from_max, ...
        NRangeAdequate, ...
        NProfileStable, ...
        ProfileJumpCountBefore, ...
        ProfileJumpCountAfter, ...
        delta_t_s, ...
        R2, ...
        RMSE, ...
        MAE, ...
        Correlation, ...
        SSE, ...
        Success, ...
        Converged, ...
        ExitFlag, ...
        Iterations, ...
        FunctionEvaluations, ...
        TotalStarts, ...
        ValidStarts, ...
        ConvergedStarts, ...
        Status);
end

%% =========================================================
%  Local function: cue-level parameter summary table
%% =========================================================

function summaryTable = create_parameter_summary_table( ...
    results, parameterNames)

    directionNames = {'Right', 'Left'};
    allNames = [parameterNames, {'N_coup'}];

    numberOfRows = ...
        numel(directionNames) * numel(allNames);

    Direction = strings(numberOfRows, 1);
    Parameter = strings(numberOfRows, 1);
    NumberOfValidTrials = zeros(numberOfRows, 1);
    Median = NaN(numberOfRows, 1);
    Q1 = NaN(numberOfRows, 1);
    Q3 = NaN(numberOfRows, 1);
    IQR = NaN(numberOfRows, 1);
    Minimum = NaN(numberOfRows, 1);
    Maximum = NaN(numberOfRows, 1);

    rowIndex = 0;

    for d = 1:numel(directionNames)

        directionName = directionNames{d};

        groupMask = ...
            [results.success].' & ...
            strcmp({results.direction}.', directionName);

        groupResults = results(groupMask);

        for p = 1:numel(allNames)

            rowIndex = rowIndex + 1;
            Direction(rowIndex) = string(directionName);
            Parameter(rowIndex) = string(allNames{p});

            if isempty(groupResults)
                continue
            end

            if p <= numel(parameterNames)

                parameterMatrix = ...
                    vertcat(groupResults.parameters);
                values = parameterMatrix(:, p);

            else

                values = [groupResults.N_coup].';
            end

            validMask = isfinite(values);
            values = values(validMask);

            NumberOfValidTrials(rowIndex) = numel(values);

            if isempty(values)
                continue
            end

            [Median(rowIndex), Q1(rowIndex), Q3(rowIndex), ...
                Minimum(rowIndex), Maximum(rowIndex)] = ...
                robust_summary(values);

            IQR(rowIndex) = Q3(rowIndex) - Q1(rowIndex);
        end
    end

    summaryTable = table( ...
        Direction, ...
        Parameter, ...
        NumberOfValidTrials, ...
        Median, ...
        Q1, ...
        Q3, ...
        IQR, ...
        Minimum, ...
        Maximum);
end

%% =========================================================
%  Local function: plot one cue direction
%% =========================================================

function figureHandle = plot_direction_fits( ...
    results, directionName, figureNumber)

    directionMask = strcmp({results.direction}, directionName);
    groupResults = results(directionMask);

    numberOfTrials = numel(groupResults);
    numberOfColumns = 4;
    numberOfRows = ceil(numberOfTrials / numberOfColumns);

    figureHandle = prepare_figure( ...
        figureNumber, ...
        sprintf('%s-cue single-trial fits', directionName));

    fitLayout = tiledlayout( ...
        numberOfRows, ...
        numberOfColumns, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    for i = 1:numberOfTrials

        nexttile(fitLayout);
        hold on;
        grid on;
        box on;

        oneResult = groupResults(i);

        if oneResult.success

            plot( ...
                oneResult.t_data, ...
                oneResult.y_data, ...
                'k.', ...
                'MarkerSize', 9);

            plot( ...
                oneResult.t_data, ...
                oneResult.y_model, ...
                'r-', ...
                'LineWidth', 1.5);

            xline( ...
                oneResult.parameters(7), ...
                'b--', ...
                'LineWidth', 0.8);

            title(sprintf( ...
                'Trial %d, R^2 = %.4f', ...
                oneResult.trial_number, ...
                oneResult.metrics.R2), 'Interpreter', 'tex');

        else

            text( ...
                0.5, ...
                0.5, ...
                sprintf('Trial %d\\nFit failed', ...
                    oneResult.trial_number), ...
                'Units', 'normalized', ...
                'HorizontalAlignment', 'center');

            title(sprintf('Trial %d', oneResult.trial_number));
        end

        ylim([-0.05, 1.10]);
    end

    % Shared labels appear once per figure, outside the individual tiles.
    xlabel(fitLayout, 'Cue-locked time (s)');
    ylabel(fitLayout, 'Normalized IMF5');

    sgtitle(fitLayout, sprintf( ...
        'Data A %s-cue single-trial original-model fits', ...
        directionName));
end

%% =========================================================
%  Local function: plot parameter distributions
%% =========================================================

function figureHandle = plot_parameter_distributions( ...
    results, parameterNames, figureNumber)

    allNames = [parameterNames, {'N_coup'}];

    figureHandle = prepare_figure( ...
        figureNumber, ...
        'Parameter distributions');

    tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

    for p = 1:numel(allNames)

        [rightValues, leftValues] = ...
            get_parameter_values(results, p);

        nexttile;
        plot_two_group_distribution( ...
            rightValues, ...
            leftValues, ...
            allNames{p}, ...
            allNames{p});
    end

    sgtitle( ...
        'Data A parameter distributions: points, median, and IQR');
end

%% =========================================================
%  Local function: create a numbered figure compatibly
%% =========================================================

function figureHandle = prepare_figure(figureNumber, figureName)

    % MATLAB does not accept figure(number, 'Property', value).
    % Create/select the numbered figure first, then set its properties.
    figureHandle = figure(figureNumber);
    clf(figureHandle);

    set( ...
        figureHandle, ...
        'Color', 'white', ...
        'Name', figureName, ...
        'NumberTitle', 'off');
end

%% =========================================================
%  Local function: retrieve parameter values for plotting
%% =========================================================

function [rightValues, leftValues] = ...
    get_parameter_values(results, parameterIndex)

    rightValues = [];
    leftValues = [];

    for i = 1:numel(results)

        if ~results(i).success
            continue
        end

        if parameterIndex <= numel(results(i).parameters)
            value = results(i).parameters(parameterIndex);
        else
            value = results(i).N_coup;
        end

        if strcmp(results(i).direction, 'Right')
            rightValues(end + 1, 1) = value; %#ok<AGROW>
        else
            leftValues(end + 1, 1) = value; %#ok<AGROW>
        end
    end
end

%% =========================================================
%  Local functions: compact distribution plot
%% =========================================================

function plot_two_group_distribution( ...
    rightValues, leftValues, yLabelText, titleText)

    hold on;
    grid on;
    box on;

    plot_one_distribution( ...
        1, rightValues, [0.85, 0.25, 0.20]);

    plot_one_distribution( ...
        2, leftValues, [0.20, 0.40, 0.85]);

    xlim([0.5, 2.5]);
    set(gca, ...
        'XTick', 1:2, ...
        'XTickLabel', {'Right', 'Left'});

    ylabel(yLabelText, 'Interpreter', 'none');
    title(titleText, 'Interpreter', 'none');
end

function plot_one_distribution(xPosition, values, plotColor)

    values = values(:);
    values = sort(values(isfinite(values)));

    if isempty(values)
        return
    end

    if numel(values) == 1
        jitter = 0;
    else
        jitter = linspace(-0.08, 0.08, numel(values)).';
    end

    scatter( ...
        xPosition + jitter, ...
        values, ...
        30, ...
        plotColor, ...
        'filled');

    [medianValue, q1Value, q3Value] = ...
        robust_summary(values);

    plot( ...
        [xPosition, xPosition], ...
        [q1Value, q3Value], ...
        '-', ...
        'Color', plotColor, ...
        'LineWidth', 6);

    plot( ...
        [xPosition - 0.13, xPosition + 0.13], ...
        [medianValue, medianValue], ...
        'k-', ...
        'LineWidth', 2);
end

%% =========================================================
%  Local functions: robust summaries
%% =========================================================

function [medianValue, q1Value, q3Value, ...
    minimumValue, maximumValue] = robust_summary(values)

    values = sort(values(isfinite(values)));
    values = values(:);

    if isempty(values)
        medianValue = NaN;
        q1Value = NaN;
        q3Value = NaN;
        minimumValue = NaN;
        maximumValue = NaN;
        return
    end

    medianValue = linear_percentile(values, 50);
    q1Value = linear_percentile(values, 25);
    q3Value = linear_percentile(values, 75);
    minimumValue = values(1);
    maximumValue = values(end);
end

function value = linear_percentile(sortedValues, percentage)

    sortedValues = sort(sortedValues(:));
    numberOfValues = numel(sortedValues);

    if numberOfValues == 1
        value = sortedValues(1);
        return
    end

    position = ...
        1 + (numberOfValues - 1) * percentage / 100;

    lowerIndex = floor(position);
    upperIndex = ceil(position);

    if lowerIndex == upperIndex
        value = sortedValues(lowerIndex);
    else
        weight = position - lowerIndex;
        value = ...
            (1 - weight) * sortedValues(lowerIndex) + ...
            weight * sortedValues(upperIndex);
    end
end

%% =========================================================
%  Local functions: progress and utilities
%% =========================================================

function print_trial_progress( ...
    stageName, currentIndex, totalCount, result, ...
    elapsedTrial, remainingSeconds)

    if result.success
        if result.N_range_adequate
            nRangeText = 'N range OK';
        else
            nRangeText = 'CHECK N upper';
        end

        if isfinite(result.profile_jump_count_after)
            profileText = sprintf( ...
                'jumps %d->%d', ...
                result.profile_jump_count_before, ...
                result.profile_jump_count_after);
        else
            profileText = 'coarse profile';
        end

        fprintf([ ...
            '%s %2d/%2d | %s trial %d | N %d | ', ...
            'R^2 %.4f | RMSE %.4f | MAE %.4f | r %.4f | ', ...
            'SSE %.4g | converged %d | ExitFlag %g | %s | %s | ', ...
            '%.1f s | ETA %.1f min\n'], ...
            stageName, ...
            currentIndex, ...
            totalCount, ...
            result.direction, ...
            result.trial_number, ...
            result.N_coup, ...
            result.metrics.R2, ...
            result.metrics.RMSE, ...
            result.metrics.MAE, ...
            result.metrics.Correlation, ...
            result.metrics.SSE, ...
            result.converged, ...
            result.exitflag, ...
            nRangeText, ...
            profileText, ...
            elapsedTrial, ...
            remainingSeconds / 60);
    else
        fprintf([ ...
            '%s %2d/%2d | %s trial %d | FAILED | %s\n'], ...
            stageName, ...
            currentIndex, ...
            totalCount, ...
            result.direction, ...
            result.trial_number, ...
            result.status_message);
    end
end

function value = get_output_value(outputStructure, fieldName, defaultValue)

    if isstruct(outputStructure) && ...
       isfield(outputStructure, fieldName)
        value = outputStructure.(fieldName);
    else
        value = defaultValue;
    end
end

function value = get_optional_field(structure, fieldName, defaultValue)

    if isfield(structure, fieldName)
        value = structure.(fieldName);
    else
        value = defaultValue;
    end
end

%% =========================================================
%  Local function: empty metric structure
%% =========================================================

function metrics = empty_metrics()

    metrics = struct;
    metrics.R2 = -inf;
    metrics.RMSE = inf;
    metrics.MAE = inf;
    metrics.Correlation = NaN;
    metrics.SSE = inf;
end

%% =========================================================
%  Local function: empty per-N candidate
%% =========================================================

function candidate = empty_candidate()

    candidate = struct;
    candidate.success = false;
    candidate.N_coup = NaN;
    candidate.parameters = NaN(1, 7);
    candidate.y_model = [];
    candidate.metrics = empty_metrics();
    candidate.R2 = -inf;
    candidate.RMSE = inf;
    candidate.exitflag = NaN;
    candidate.converged = false;
    candidate.iterations = NaN;
    candidate.function_evaluations = NaN;
    candidate.first_order_optimality = NaN;
    candidate.total_starts = 0;
    candidate.valid_starts = 0;
    candidate.converged_starts = 0;
    candidate.rescue_count = 0;
end

%% =========================================================
%  Local function: empty trial result
%% =========================================================

function result = empty_trial_result()

    result = struct;
    result.direction = '';
    result.trial_number = NaN;
    result.cue_sample = NaN;
    result.raw_peak_amplitude = NaN;
    result.peak_time_cue_locked = NaN;
    result.fit_start_time_cue_locked = NaN;
    result.fit_end_time_cue_locked = NaN;

    result.t_data = [];
    result.y_data = [];
    result.number_of_samples = 0;

    result.success = false;
    result.status_message = 'Not fitted';
    result.parameters = NaN(1, 7);
    result.N_coup = NaN;
    result.y_model = [];
    result.residual = [];
    result.metrics = empty_metrics();

    result.exitflag = NaN;
    result.converged = false;
    result.iterations = NaN;
    result.function_evaluations = NaN;
    result.first_order_optimality = NaN;

    result.total_starts = 0;
    result.valid_starts = 0;
    result.converged_starts = 0;

    result.lower_bounds = NaN(1, 7);
    result.upper_bounds = NaN(1, 7);

    result.N_at_lower_candidate = false;
    result.N_at_upper_candidate = false;
    result.N_of_max_R2 = NaN;
    result.max_R2_across_N = NaN;
    result.R2_loss_from_max = NaN;
    result.N_range_adequate = false;
    result.plateau_R2_tolerance = NaN;
    result.plateau_reference_R2 = NaN;
    result.N_profile_stable = false;
    result.profile_jump_count_before = NaN;
    result.profile_jump_count_after = NaN;
    result.profile_drop_tolerance = NaN;
    result.per_N = repmat(empty_candidate(), 0, 1);
end
