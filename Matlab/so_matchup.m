function [T_train, T_val] = so_matchup(csvFile, relErrMax, trainFrac)
% SO_MATCHUP
%   Load matchup CSV, filter by relative error between OCChla and Tchla,
%   and perform a year-stratified train/validation split.
%
% Inputs:
%   csvFile   - path to 'matchup_dataset_SO_zenodo.csv'
%   relErrMax - maximum relative error (%) |OCChla - Tchla|/Tchla*100 (default 50)
%   trainFrac - fraction per year used for training (default 0.70)
%
% Outputs:
%   T_train   - table of training samples
%   T_val     - table of validation samples
%
% Notes:
%   - No intermediate CSV files are written.
%   - All downstream analysis should work directly from T_train and T_val.

    if nargin < 2 || isempty(relErrMax)
        relErrMax = 50;       % default: allow up to 50% mismatch
    end
    if nargin < 3 || isempty(trainFrac)
        trainFrac = 0.70;     % default: 70% training, 30% validation per year
    end

    % -------------------- Load full matchup dataset --------------------
    T = readtable(csvFile);

    % Basic checks
    requiredVars = {'Tchla','OCChla','Year'};
    if ~all(ismember(requiredVars, T.Properties.VariableNames))
        error('Input table must contain columns: %s', strjoin(requiredVars, ', '));
    end

    Tchla  = T.Tchla;
    OCChla = T.OCChla;

    % -------------------- Relative error filter ------------------------
    maskValid = isfinite(Tchla) & isfinite(OCChla) & Tchla > 0;
    percentErr = NaN(size(Tchla));
    percentErr(maskValid) = abs(OCChla(maskValid) - Tchla(maskValid)) ./ Tchla(maskValid) * 100;

    keep = maskValid & (percentErr <= relErrMax);
    T = T(keep, :);

    fprintf('After %.1f%% error filter: %d matchup samples retained.\n', relErrMax, height(T));

    % -------------------- Year-stratified split ------------------------
    years = unique(T.Year);
    T_train = T([],:);   % empty tables with same variables
    T_val   = T([],:);

    rng(1);              % reproducible

    for i = 1:numel(years)
        yr = years(i);
        idxYear = find(T.Year == yr);
        N = numel(idxYear);

        if N <= 1
            % too few samples to split; keep all in training
            T_train = [T_train; T(idxYear,:)];
            continue;
        end

        idxRand = idxYear(randperm(N));
        N_train = max(1, round(trainFrac * N));

        trainIdx = idxRand(1:N_train);
        valIdx   = idxRand(N_train+1:end);

        T_train = [T_train; T(trainIdx,:)];
        T_val   = [T_val;   T(valIdx,:)];
    end

    fprintf('Training set size : %d rows\n', height(T_train));
    fprintf('Validation set size: %d rows\n', height(T_val));
end
