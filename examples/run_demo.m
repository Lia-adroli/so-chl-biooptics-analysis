%% ============================================================
%  RUN_DEMO — Southern Ocean Bio-Optical Model (Demo Example)
%
%  This script demonstrates the full PSC/PFT reconstruction
%  workflow using a small example matchup dataset.
%
%  It performs:
%   1) Data filtering and year-stratified split
%   2) NNLS + GBM training with K-fold cross-validation
%   3) Independent validation using the bin-based model
%
%  This demo is for testing and illustration only.
%  Scientific results use the full Zenodo dataset:
%  DOI: 10.5281/zenodo.17875100
% ============================================================

%% --- Setup paths ---
addpath('../Matlab')
addpath('../examples')

%% --- Input demo data ---
csvFile   = 'demo_matchup.csv';   % demo input file
relErrMax = 50;                   % maximum OCChla–Tchla error (%)
trainFrac = 0.70;                 % 70% training, 30% validation
Kfold     = 3;                    % K-fold cross-validation

%% --- Step 1: Prepare matchups (filter + split) ---
[T_train, T_val] = so_matchup(csvFile, relErrMax, trainFrac);

%% --- Step 2: Train model (NNLS + GBM + bin model) ---
model = train_model(T_train, Kfold);

%% --- Step 3: Independent validation ---
stats = validate_model(T_val, model);

%% --- Display validation results ---
disp(' ')
disp('===== VALIDATION STATISTICS (LOG10 SPACE) =====')
for i = 1:numel(stats)
    fprintf('%-20s  N=%4d  r=%.3f  R2=%.3f  RMSE=%.3f  MAE=%.3f\n', ...
        stats(i).name, ...
        stats(i).n, ...
        stats(i).r, ...
        stats(i).R2, ...
        stats(i).RMSE, ...
        stats(i).MAE);
end

disp(' ')
disp('Demo run completed successfully.')
