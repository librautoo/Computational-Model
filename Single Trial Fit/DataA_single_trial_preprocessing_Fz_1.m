%% Data A: single-trial preprocessing of Fz IMF5
% This script:
%   1. uses only the manually retained left- and right-cue trials;
%   2. extracts the 400 samples immediately after each cue;
%   3. applies VMD independently to every Fz trial;
%   4. extracts IMF5 independently from every trial;
%   5. normalizes each IMF5 by its own maximum positive peak;
%   6. saves the original IMF5 peak amplitude before normalization;
%   7. extracts the continuous positive interval containing that peak;
%   8. preserves the original cue-locked time for later model fitting;
%   9. creates separate right- and left-cue IMF5 overview figures;
%  10. creates one combined cue-locked positive fitting-window figure.
%
% Important:
% The fitting time is NOT shifted so that the positive-wave onset becomes
% t = 0. All saved times remain relative to the original cue onset.

clearvars;
close all;
clc;

%% =========================================================
%  File names
%% =========================================================

inputFile = 'yu20150723a_Task-1.mat';
outputMatFile = 'DataA_Fz_single_trial_preprocessed.mat';
outputCsvFile = 'DataA_Fz_single_trial_peak_summary.csv';

%% =========================================================
%  Fixed preprocessing settings
%% =========================================================

fs = 256;

FzRow = 1;
CueMarkerRow = 8;

NumIMFs = 6;
PenaltyFactor = 2000;
MaxIterations = 500;
TargetIMF = 5;

Len = 400;

% This range is used only as a quality-control flag.
% It does not determine the peak and does not remove any trial.
ExpectedPeakRange = [0.5, 0.9];

% The first extracted point is one sample after cue onset.
% Therefore, it occurs at 1/fs seconds rather than exactly 0 seconds.
cueLockedTime = (1:Len) / fs;

%% =========================================================
%  Manually retained Data A trials
%% =========================================================

keepLeft = [9, 14, 19, 27, 35, 37, 38, 39, 40, 41, 42, 46, 47, 54];
keepRight = [8, 17, 18, 22, 23, 24, 29, 33, 37, 39, 40, 42, 43, 44];

%% =========================================================
%  Load Data A
%% =========================================================

if ~isfile(inputFile)
    error('Input file not found: %s', inputFile);
end

loadedData = load(inputFile);

if ~isfield(loadedData, 'data')
    error('The input MAT file does not contain a variable named "data".');
end

data = loadedData.data;
[numberOfRows, numberOfSamples] = size(data);

if numberOfRows < max(FzRow, CueMarkerRow)
    error('The data matrix does not contain the required Fz and cue-marker rows.');
end

%% =========================================================
%  Detect cue onsets
%% =========================================================

cueSignal = data(CueMarkerRow, :);

LeftwardsCue = find( ...
    cueSignal(2:end) == -1 & cueSignal(1:end-1) == 0) + 1;

RightwardsCue = find( ...
    cueSignal(2:end) == 1 & cueSignal(1:end-1) == 0) + 1;

numberOfLeftCues = numel(LeftwardsCue);
numberOfRightCues = numel(RightwardsCue);

fprintf('Task-1 Data A: total left-cue trials  = %d\n', numberOfLeftCues);
fprintf('Task-1 Data A: total right-cue trials = %d\n', numberOfRightCues);

if any(keepLeft < 1) || any(keepLeft > numberOfLeftCues)
    error('At least one value in keepLeft is outside the available left-trial range.');
end

if any(keepRight < 1) || any(keepRight > numberOfRightCues)
    error('At least one value in keepRight is outside the available right-trial range.');
end

%% =========================================================
%  Store settings
%% =========================================================

settings = struct;
settings.fs = fs;
settings.FzRow = FzRow;
settings.CueMarkerRow = CueMarkerRow;
settings.NumIMFs = NumIMFs;
settings.PenaltyFactor = PenaltyFactor;
settings.MaxIterations = MaxIterations;
settings.TargetIMF = TargetIMF;
settings.WindowLengthSamples = Len;
settings.WindowDurationSeconds = Len / fs;
settings.ExpectedPeakRangeSeconds = ExpectedPeakRange;
settings.FirstExtractedSampleAfterCue = 1;

%% =========================================================
%  Process every retained trial independently
%% =========================================================

rightTrials = preprocess_one_direction( ...
    data, RightwardsCue, keepRight, 1, 'Right', ...
    cueLockedTime, settings, numberOfSamples);

leftTrials = preprocess_one_direction( ...
    data, LeftwardsCue, keepLeft, -1, 'Left', ...
    cueLockedTime, settings, numberOfSamples);

%% =========================================================
%  Create a compact summary table
%% =========================================================

rightSummary = create_summary_table(rightTrials);
leftSummary = create_summary_table(leftTrials);
summaryTable = [rightSummary; leftSummary];

disp(' ');
disp('Single-trial preprocessing summary:');
disp(summaryTable);

%% =========================================================
%  Collect and save all outputs
%% =========================================================

DataA = struct;
DataA.dataset = 'yu20150723a_Task-1';
DataA.channel = 'Fz';
DataA.description = ...
    ['Single-trial IMF5 preprocessing with original cue-locked time. ', ...
     'No across-trial averaging and no time-origin shifting were applied.'];
DataA.settings = settings;
DataA.cue_locked_time_full = cueLockedTime;
DataA.keepLeft = keepLeft;
DataA.keepRight = keepRight;
DataA.allLeftCueSamples = LeftwardsCue;
DataA.allRightCueSamples = RightwardsCue;
DataA.right = rightTrials;
DataA.left = leftTrials;
DataA.summary = summaryTable;

save(outputMatFile, 'DataA', '-v7.3');
writetable(summaryTable, outputCsvFile);

fprintf('\nSaved MATLAB results to: %s\n', outputMatFile);
fprintf('Saved peak summary to: %s\n', outputCsvFile);

%% =========================================================
%  Create quality-control figures
%% =========================================================

% The original four-panel Figure 1 is separated into two two-panel figures.
% This makes the right- and left-cue results easier to insert into LaTeX.
plot_group_imf5_overview( ...
    rightTrials, cueLockedTime, ExpectedPeakRange, 1);

plot_group_imf5_overview( ...
    leftTrials, cueLockedTime, ExpectedPeakRange, 2);

% Keep the original fitting-window summary unchanged in content.
plot_positive_fitting_windows(rightTrials, leftTrials, 3);

%% =========================================================
%  Local function: preprocess one cue direction
%% =========================================================

function group = preprocess_one_direction( ...
    data, allCueSamples, retainedTrialNumbers, cueValue, directionLabel, ...
    cueLockedTime, settings, numberOfSamples)

    numberOfTrials = numel(retainedTrialNumbers);
    Len = settings.WindowLengthSamples;

    group = struct;
    group.direction = directionLabel;
    group.cue_value = cueValue;
    group.trial_number = retainedTrialNumbers(:);
    group.cue_sample = zeros(numberOfTrials, 1);
    group.window_start_sample = zeros(numberOfTrials, 1);
    group.window_end_sample = zeros(numberOfTrials, 1);

    group.raw_fz = zeros(numberOfTrials, Len);
    group.imf5_raw = zeros(numberOfTrials, Len);
    group.imf5_normalized = zeros(numberOfTrials, Len);

    group.peak_amplitude_raw_imf5 = zeros(numberOfTrials, 1);
    group.peak_index_in_window = zeros(numberOfTrials, 1);
    group.peak_time_cue_locked = zeros(numberOfTrials, 1);
    group.peak_in_expected_range = false(numberOfTrials, 1);

    group.fit_indices = cell(numberOfTrials, 1);
    group.fit_time_cue_locked = cell(numberOfTrials, 1);
    group.fit_data_raw_imf5 = cell(numberOfTrials, 1);
    group.fit_data_normalized = cell(numberOfTrials, 1);
    group.fit_start_time_cue_locked = zeros(numberOfTrials, 1);
    group.fit_end_time_cue_locked = zeros(numberOfTrials, 1);
    group.fit_sample_count = zeros(numberOfTrials, 1);

    for i = 1:numberOfTrials

        trialNumber = retainedTrialNumbers(i);
        cueSample = allCueSamples(trialNumber);

        windowStart = cueSample + 1;
        windowEnd = cueSample + Len;

        if windowEnd > numberOfSamples
            error(['%s trial %d does not contain %d complete samples ', ...
                   'after cue onset.'], ...
                  directionLabel, trialNumber, Len);
        end

        rawFz = data(settings.FzRow, windowStart:windowEnd);

        [imf, ~, ~] = vmd_one_channel( ...
            rawFz, ...
            settings.NumIMFs, ...
            settings.PenaltyFactor, ...
            settings.MaxIterations);

        if size(imf, 2) < settings.TargetIMF
            error('%s trial %d returned fewer than %d IMFs.', ...
                  directionLabel, trialNumber, settings.TargetIMF);
        end

        imf5 = imf(:, settings.TargetIMF).';

        % Use the maximum positive IMF5 peak, as in the previous analysis.
        [peakAmplitude, peakIndex] = max(imf5);

        if ~isfinite(peakAmplitude) || peakAmplitude <= 0
            error(['%s trial %d has no valid positive IMF5 peak and ', ...
                   'cannot be normalized by a positive peak.'], ...
                  directionLabel, trialNumber);
        end

        normalizedImf5 = imf5 / peakAmplitude;

        % Find the full continuous positive interval containing the peak.
        % The search uses the complete 400-sample IMF5 waveform.
        fitStartIndex = peakIndex;
        while fitStartIndex > 1 && normalizedImf5(fitStartIndex - 1) > 0
            fitStartIndex = fitStartIndex - 1;
        end

        fitEndIndex = peakIndex;
        while fitEndIndex < Len && normalizedImf5(fitEndIndex + 1) > 0
            fitEndIndex = fitEndIndex + 1;
        end

        fitIndices = fitStartIndex:fitEndIndex;

        group.cue_sample(i) = cueSample;
        group.window_start_sample(i) = windowStart;
        group.window_end_sample(i) = windowEnd;

        group.raw_fz(i, :) = rawFz;
        group.imf5_raw(i, :) = imf5;
        group.imf5_normalized(i, :) = normalizedImf5;

        group.peak_amplitude_raw_imf5(i) = peakAmplitude;
        group.peak_index_in_window(i) = peakIndex;
        group.peak_time_cue_locked(i) = cueLockedTime(peakIndex);
        group.peak_in_expected_range(i) = ...
            cueLockedTime(peakIndex) >= settings.ExpectedPeakRangeSeconds(1) && ...
            cueLockedTime(peakIndex) <= settings.ExpectedPeakRangeSeconds(2);

        group.fit_indices{i} = fitIndices;

        % Keep the original cue-locked time.
        % Do not subtract cueLockedTime(fitStartIndex).
        group.fit_time_cue_locked{i} = cueLockedTime(fitIndices);

        group.fit_data_raw_imf5{i} = imf5(fitIndices);
        group.fit_data_normalized{i} = normalizedImf5(fitIndices);
        group.fit_start_time_cue_locked(i) = cueLockedTime(fitStartIndex);
        group.fit_end_time_cue_locked(i) = cueLockedTime(fitEndIndex);
        group.fit_sample_count(i) = numel(fitIndices);

        fprintf(['%s trial %2d | cue sample = %7d | raw IMF5 peak = % .6g ', ...
                 '| peak time = %.4f s | positive window = %.4f to %.4f s ', ...
                 '| samples = %d\n'], ...
                directionLabel, ...
                trialNumber, ...
                cueSample, ...
                peakAmplitude, ...
                cueLockedTime(peakIndex), ...
                cueLockedTime(fitStartIndex), ...
                cueLockedTime(fitEndIndex), ...
                numel(fitIndices));
    end
end

%% =========================================================
%  Local function: VMD of one channel and one trial
%% =========================================================

function [imf, residual, info] = ...
    vmd_one_channel(x, NumIMFs, PenaltyFactor, MaxIterations)

    x = x(:);

    [imf, residual, info] = vmd(x, ...
        NumIMFs = NumIMFs, ...
        PenaltyFactor = PenaltyFactor, ...
        MaxIterations = MaxIterations);

    % Ensure that rows are time points and columns are IMFs.
    if size(imf, 1) ~= length(x) && size(imf, 2) == length(x)
        imf = imf.';
    end

    if size(imf, 1) ~= length(x)
        error('Unexpected IMF array size returned by vmd.');
    end
end

%% =========================================================
%  Local function: build the summary table
%% =========================================================

function summaryTable = create_summary_table(group)

    numberOfTrials = numel(group.trial_number);
    directionColumn = repmat(string(group.direction), numberOfTrials, 1);

    summaryTable = table( ...
        directionColumn, ...
        group.trial_number, ...
        group.cue_sample, ...
        group.window_start_sample, ...
        group.window_end_sample, ...
        group.peak_amplitude_raw_imf5, ...
        group.peak_time_cue_locked, ...
        group.peak_in_expected_range, ...
        group.fit_start_time_cue_locked, ...
        group.fit_end_time_cue_locked, ...
        group.fit_sample_count, ...
        'VariableNames', { ...
            'Direction', ...
            'TrialNumber', ...
            'CueSample', ...
            'WindowStartSample', ...
            'WindowEndSample', ...
            'RawPeakAmplitudeIMF5', ...
            'PeakTimeCueLocked_s', ...
            'PeakInExpectedRange', ...
            'FitStartCueLocked_s', ...
            'FitEndCueLocked_s', ...
            'FitSampleCount'});
end

%% =========================================================
%  Local function: plot one cue direction as a two-panel figure
%% =========================================================

function plot_group_imf5_overview( ...
    group, cueLockedTime, ExpectedPeakRange, figureNumber)

    numberOfTrials = numel(group.trial_number);
    trialColors = lines(numberOfTrials);

    figure(figureNumber);
    clf;
    set(gcf, ...
        'Color', 'w', ...
        'Name', sprintf('Data A Fz - %s cue IMF5 overview', group.direction), ...
        'NumberTitle', 'off');

    tiledlayout(2, 1, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    % Upper panel: original IMF5 waveforms and their positive peaks.
    nexttile;
    hold on;

    for i = 1:numberOfTrials
        plot( ...
            cueLockedTime, ...
            group.imf5_raw(i, :), ...
            'Color', trialColors(i, :), ...
            'LineWidth', 1.0);

        scatter( ...
            group.peak_time_cue_locked(i), ...
            group.peak_amplitude_raw_imf5(i), ...
            24, ...
            trialColors(i, :), ...
            'filled');
    end

    yline(0, 'k:');
    xline(ExpectedPeakRange(1), 'k--');
    xline(ExpectedPeakRange(2), 'k--');
    xlabel('Time after cue onset (s)');
    ylabel('Original IMF5 amplitude');
    title(sprintf( ...
        '%s cue: original IMF5 and detected positive peaks', ...
        group.direction));
    grid on;
    box on;
    xlim([cueLockedTime(1), cueLockedTime(end)]);

    % Lower panel: independently normalized IMF5 waveforms.
    nexttile;
    hold on;

    for i = 1:numberOfTrials
        plot( ...
            cueLockedTime, ...
            group.imf5_normalized(i, :), ...
            'Color', trialColors(i, :), ...
            'LineWidth', 1.0);

        scatter( ...
            group.peak_time_cue_locked(i), ...
            1, ...
            24, ...
            trialColors(i, :), ...
            'filled');
    end

    yline(0, 'k:');
    xline(ExpectedPeakRange(1), 'k--');
    xline(ExpectedPeakRange(2), 'k--');
    xlabel('Time after cue onset (s)');
    ylabel('Normalized IMF5 amplitude');
    title(sprintf( ...
        '%s cue: IMF5 normalized by each trial''s own positive peak', ...
        group.direction));
    grid on;
    box on;
    xlim([cueLockedTime(1), cueLockedTime(end)]);

    sgtitle(sprintf( ...
        'Data A Fz single-trial IMF5 quality control: %s cue (%d trials)', ...
        group.direction, ...
        numberOfTrials));
end

%% =========================================================
%  Local function: plot all cue-locked positive fitting windows
%% =========================================================

function plot_positive_fitting_windows( ...
    rightGroup, leftGroup, figureNumber)

    figure(figureNumber);
    clf;
    set(gcf, ...
        'Color', 'w', ...
        'Name', 'Data A cue-locked positive fitting windows', ...
        'NumberTitle', 'off');

    tiledlayout(1, 2, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    plot_one_fitting_window_group(rightGroup);
    plot_one_fitting_window_group(leftGroup);

    sgtitle('Data A cue-locked positive fitting windows');
end

%% =========================================================
%  Local function: plot one direction in the fitting-window figure
%% =========================================================

function plot_one_fitting_window_group(group)

    numberOfTrials = numel(group.trial_number);
    trialColors = lines(numberOfTrials);

    nexttile;
    hold on;

    for i = 1:numberOfTrials
        plot( ...
            group.fit_time_cue_locked{i}, ...
            group.fit_data_normalized{i}, ...
            'Color', trialColors(i, :), ...
            'LineWidth', 1.2);
    end

    xlabel('Time after cue onset (s)');
    ylabel('Normalized amplitude');
    title(sprintf( ...
        '%s cue: %d retained trials', ...
        group.direction, ...
        numberOfTrials));
    grid on;
    box on;
end
