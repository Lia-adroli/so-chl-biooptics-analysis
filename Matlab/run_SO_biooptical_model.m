function results = run_SO_biooptical_model(csvFile, varargin)
% RUN_SO_BIOOPTICAL_MODEL
%   End-to-end Southern Ocean bio-optical PSC/PFT workflow.
%
%   This function:
%     1) Loads the matchup dataset (HPLC + OC-CCI)
%     2) Applies relative-error filtering between OCChla and Tchla
%     3) Performs a year-stratified train/validation split
%     4) Trains the NNLS + GBM + K-fold cross-validated model
%     5) Builds a bin-based empirical PSC/PFT model
%     6) Validates against independent OC-CCI–based reconstruction
%
% SYNTAX
%   results = run_SO_biooptical_model(csvFile)
%   results = run_SO_biooptical_model(csvFile, 'RelErrMax', 50, ...
%                                      'TrainFrac', 0.7, ...
%                                      'KFolds', 5)
%
% INPUTS
%   csvFile   - path to 'matchup_dataset_SO_zenodo.csv'
%
% NAME–VALUE PAIRS (optional)
%   'RelErrMax'  - maximum relative error (%) between OCChla and Tchla
%                  used in initial filtering (default: 50)
%   'TrainFrac'  - fraction (0–1) of data per year used for training
%                  in year-stratified split (default: 0.70)
%   'KFolds'     - number of folds for GBM K-fold cross-validation
%                  (default: 5)
%
% OUTPUT
%   results - struct with fields:
%       .T_train  - training table after filtering & split
%       .T_val    - validation table after filtering & split
%       .model    - trained model struct from train_model.m
%       .stats    - validation skill metrics from validate_model.m
%
% EXAMPLE
%   results = run_SO_biooptical_model('data/matchup_dataset_SO_zenodo.csv');
%
%   % Access bin model coefficients:
%   coeffs = results.model.binModel.coeffs;
%
%   % Access validation R² for PSC_micro:
%   idx = strcmp({results.stats.name}, 'PSC_micro');
%   results.stats(idx).R2
%
% REQUIREMENTS
%   The following functions must be on the MATLAB path:
%     - so_matchups.m
%     - train_model.m
%     - bin_model.m
%     - validate_model.m
%
%   The input CSV must contain (at minimum) these columns:
%     Tchla, OCChla, Year, lat,
%     Fuco, Per, X19hex, Allo, X19but, Chl_b, Zea, DVChla
%
% AUTHOR
%   (c) Your Name / Lia Adroli, YEAR
%   Southern Ocean PSC/PFT bio-optical reconstruction framework.

    %% -------------------- Defaults & parsing ---------------------------
    if nargin < 1 || isempty(csvFile)
        error('csvFile (path to matchup_dataset_SO_zenodo.csv) must be provided.');
    end

    % Default options
    opts.RelErrMax = 50;   % % relative error filter between OCChla & Tchla
    opts.TrainFrac = 0.70; % year-wise training fraction
    opts.KFolds    = 5;    % K-fold cross-validation

    % Parse name–value pairs
    if mod(numel(varargin),2) ~= 0
        error('Optional arguments must be provided as name–value pairs.');
    end

    for k = 1:2:numel(varargin)
        name  = varargin{k};
        value = varargin{k+1};

        if ~ischar(name) && ~isstring(name)
            error('Optional argument name must be a character vector or string.');
        end

        switch lower(strtrim(name))
            case 'relerrmax'
                opts.RelErrMax = value;
            case 'trainfrac'
                opts.TrainFrac = value;
            case 'kfolds'
                opts.KFolds = value;
            otherwise
                error('Unrecognized option name: %s', name);
        end
    end

    %% -------------------- 1) Prepare matchups -------------------------
    fprintf('=== STEP 1: Loading & filtering matchup dataset ===\n');
    fprintf('  File          : %s\n', csvFile);
    fprintf('  RelErrMax     : %.1f %%\n', opts.RelErrMax);
    fprintf('  TrainFrac     : %.2f (year-stratified)\n', opts.TrainFrac);

    [T_train, T_val] = so_matchups(csvFile, opts.RelErrMax, opts.TrainFrac);

    %% -------------------- 2) Train bio-optical model ------------------
    fprintf('\n=== STEP 2: Training NNLS + GBM + bin model ===\n');
    fprintf('  K-fold CV     : %d folds\n', opts.KFolds);

    model = train_model(T_train, opts.KFolds);

    %% -------------------- 3) Independent validation -------------------
    fprintf('\n=== STEP 3: Validating using OCChla (OC-CCI) ===\n');

    stats = validate_model(T_val, model);

    %% -------------------- 4) Summary to console -----------------------
    fprintf('\n=== VALIDATION SUMMARY (log10-space metrics) ===\n');
    try
        % Convert stats struct to table for a clean print
        statsTable = struct2table(stats);
        disp(statsTable);
    catch
        % Fallback: simple loop
        for i = 1:numel(stats)
            fprintf('%-18s : n = %4d | r = %.3f | R^2 = %.3f | RMSE = %.3f | MAE = %.3f\n', ...
                stats(i).name, stats(i).n, stats(i).r, stats(i).R2, ...
                stats(i).RMSE, stats(i).MAE);
        end
    end

    %% -------------------- 5) Pack outputs -----------------------------
    results = struct();
    results.T_train = T_train;
    results.T_val   = T_val;
    results.model   = model;
    results.stats   = stats;

    fprintf('\n=== DONE: run_SO_biooptical_model completed successfully ===\n');
end
