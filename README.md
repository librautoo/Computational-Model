# Computational Model

MATLAB code for the EEG preprocessing, single-trial coupling–release model fitting, and subsequent three-delay analysis used in this project.

The repository contains three main parts:

```text
Computational-Model/
├── Preprocessing/
├── Single Trial Fit/
└── 3V Delay/
```

The number at the end of each MATLAB filename (for example, `_1`, `_2`, `_3`) indicates the intended running order within that workflow.

## 1. Preprocessing

The `Preprocessing` folder contains the preprocessing and IMF5 inspection workflow for Data A and Data B.

The two original Task-1 datasets are:

```text
yu20150723a_Task-1.mat   # Data A
yu20150723b_Task-1.mat   # Data B
```

Run the corresponding scripts for Data A or Data B in numeric order.

### Data A

```text
task1_sift_a_1.m
        ↓
task1_VMD_a_2.m
        ↓
task1_IMF5_norm_a_3.m
```

### Data B

```text
task1_sift_b_1.m
        ↓
task1_VMD_b_2.m
        ↓
task1_IMF5_normalized__b_3.m
```

The VMD scripts use six modes with a penalty factor of 2000 and a maximum of 500 iterations. The analysis uses the manually retained trials and focuses on IMF5 extracted from the post-cue signal.

This folder is useful for the preprocessing/inspection workflow and for examining the retained IMF5 responses.

## 2. Single-Trial Coupling–Release Model Fit

The `Single Trial Fit` folder contains the complete single-trial analysis for the Fz channel.

Importantly, the single-trial preprocessing scripts are self-contained: **they perform VMD internally and do not require the VMD output from the separate `Preprocessing` folder.**

For every retained trial, the `_1` preprocessing script:

1. loads the original Task-1 dataset;
2. detects the cue onset;
3. extracts the first 400 post-cue samples from Fz;
4. applies VMD independently to that trial;
5. decomposes the signal into six modes;
6. extracts IMF5;
7. normalizes IMF5 by that trial's own maximum positive peak;
8. identifies the continuous positive interval containing the peak; and
9. preserves the original cue-locked time axis for subsequent model fitting.


### Data A

Run:

```text
DataA_single_trial_preprocessing_Fz_1.m
        ↓
DataA_single_trial_original_model_fit_Fz_2.m
```

The first script produces:

```text
DataA_Fz_single_trial_preprocessed.mat
DataA_Fz_single_trial_peak_summary.csv
```

The second script reads the preprocessed MAT file and performs the single-trial coupling–release model fitting. Its main MAT output is:

```text
DataA_Fz_single_trial_original_model_fit_fixed_ranges.mat
```

### Data B

Run:

```text
DataB_single_trial_preprocessing_Fz_1.m
        ↓
DataB_single_trial_original_model_fit_Fz_2.m
```

The first script produces:

```text
DataB_Fz_single_trial_preprocessed.mat
DataB_Fz_single_trial_peak_summary.csv
```

The second script produces the main model-fit MAT file:

```text
DataB_Fz_single_trial_original_model_fit_fixed_ranges.mat
```

The model-fit MAT files contain the `FitAnalysis` structure used by the three-delay analysis.

## 3. Three-Delay Analysis

After completing the single-trial model fitting for both datasets, copy these **two model-fit MAT files** into the `3V Delay` folder:

```text
DataA_Fz_single_trial_original_model_fit_fixed_ranges.mat
DataB_Fz_single_trial_original_model_fit_fixed_ranges.mat
```

The preprocessed MAT files are not the inputs to this stage; the three-delay scripts read the fitted `FitAnalysis` results.

Two complementary three-delay procedures are provided:

```text
music_auto_search_all_trials_window_valid.m
atom_splitting_auto_all_trials_noise_stability.m
```

### MUSIC

Open:

```text
music_auto_search_all_trials_window_valid.m
```

Set:

```matlab
cfg.fitFile = 'DataA_Fz_single_trial_original_model_fit_fixed_ranges.mat';
```

and run the script.

Then change the input to:

```matlab
cfg.fitFile = 'DataB_Fz_single_trial_original_model_fit_fixed_ranges.mat';
```

and run it again.

By default, the script processes both Right and Left trials and writes results under the `music_results` directory.

### Atom Splitting

Open:

```text
atom_splitting_auto_all_trials_noise_stability.m
```

First run:

```matlab
cfg.fitFile = 'DataA_Fz_single_trial_original_model_fit_fixed_ranges.mat';
```

Then change it to:

```matlab
cfg.fitFile = 'DataB_Fz_single_trial_original_model_fit_fixed_ranges.mat';
```

and run the script again.

The atom-splitting analysis uses a fixed three-atom model and includes the noise-stability analysis. Results are written under the `atom_splitting_results` directory.


## MATLAB Requirements

The code is written for MATLAB. The current scripts use functions from the following toolboxes:

- **Optimization Toolbox** — including `lsqnonlin`, `fmincon`, and `lsqnonneg`;
- **Signal Processing Toolbox** — required by the MUSIC workflow;
- MATLAB's `vmd` implementation for Variational Mode Decomposition.

The delay-analysis scripts are currently configured for the Fz channel and process both Right- and Left-cue trials by default.
