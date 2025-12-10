function model = train_model(T_train, K)
% TRAIN_MODEL
%   Train NNLS + GBM model with K-fold cross-validation and construct
%   PSC/PFT fractions and a bin-based empirical model.
%
% Inputs:
%   T_train - table from prepare_so_matchups (training subset)
%   K       - number of folds for K-fold CV (default 5)
%
% Output:
%   model   - struct containing:
%       .Ttrain        - filtered training table actually used
%       .Tchla         - in situ Tchla (vector)
%       .P             - pigment matrix [N x 8]
%       .nnlsCoeff     - NNLS pigment coefficients [8 x 1]
%       .TChla_cv      - K-fold cross-validated GBM predictions
%       .gbmFull       - GBM model trained on all filtered samples
%       .PSC           - struct of PSC fractions (0–1)
%       .PFT           - struct of PFT fractions (0–1)
%       .binModel      - struct with bin-based polynomial fits for PSC/PFT

    if nargin < 2 || isempty(K)
        K = 5;   % default K-fold
    end

    % -------------------- Basic checks --------------------
    needed = {'Tchla','lat','Fuco','Per','X19hex','Allo','X19but','Chl_b','Zea','DVChla'};
    if ~all(ismember(needed, T_train.Properties.VariableNames))
        error('T_train must contain columns: %s', strjoin(needed, ', '));
    end

    % -------------------- Extract core variables -----------------------
    Tchla = T_train.Tchla;
    lat   = T_train.lat;

    % 8 pigments
    P = [T_train.Fuco,   T_train.Per,   T_train.X19hex, ...
         T_train.Allo,   T_train.X19but, T_train.Chl_b, ...
         T_train.Zea,    T_train.DVChla];

    % -------------------- Southern Ocean latitude filter ---------------
    mask_lat = isfinite(lat) & lat >= -90 & lat <= -40;

    % -------------------- QC: Tchla & pigment availability -------------
    mask_T = isfinite(Tchla) & Tchla > 0 & Tchla < 15;

    P_meas = P;
    P_meas(~isfinite(P_meas)) = NaN;
    measured = P_meas >= 1e-3;           % pigment "measured" threshold
    has4 = sum(measured, 2) >= 4;        % at least 4 pigments measured

    keep = mask_T & has4 & mask_lat;

    % Apply mask
    T_use  = T_train(keep, :);
    Tchla  = Tchla(keep);
    P      = P(keep,:);

    % Sanitize pigment values
    P(~isfinite(P)) = 0;
    P(P < 0)        = 0;
    P(P <= 1e-3)    = 0;

    fprintf('Valid Southern Ocean samples retained for training: %d\n', numel(Tchla));

    % -------------------- NNLS pigment mixture -------------------------
    % Solve: min ||P*C - Tchla||^2, subject to C >= 0
    nnlsCoeff = lsqnonneg(P, Tchla);     % [8 x 1]
    contrib   = P .* nnlsCoeff.';       % [N x 8]
    DP_lin    = sum(contrib, 2);        % diagnostic pigment sum

    % -------------------- GBM with K-fold CV ---------------------------
    epsP  = 1e-6;
    Xfeat = [P, log(P + epsP), DP_lin, log(DP_lin + epsP)];

    tTree = templateTree('MinLeafSize',3, ...
                         'NumVariablesToSample', size(P,2));

    rng(42);                               % reproducibility
    cvp = cvpartition(numel(Tchla), 'KFold', K);

    TChla_cv = NaN(size(Tchla));

    for k = 1:K
        idxTr = training(cvp, k);
        idxTe = test(cvp, k);

        gb_k = fitrensemble(Xfeat(idxTr,:), Tchla(idxTr), ...
                            'Method','LSBoost', ...
                            'Learners', tTree, ...
                            'NumLearningCycles', 1600, ...
                            'LearnRate', 0.02, ...
                            'FResample', 0.6);

        TChla_cv(idxTe) = predict(gb_k, Xfeat(idxTe,:));
    end

    % Train a final GBM on all filtered samples (for later prediction)
    gbmFull = fitrensemble(Xfeat, Tchla, ...
                           'Method','LSBoost', ...
                           'Learners', tTree, ...
                           'NumLearningCycles', 1600, ...
                           'LearnRate', 0.02, ...
                           'FResample', 0.6);

    % -------------------- PSC/PFT fractions (GBM-anchored) -------------
    Tg = max(TChla_cv, 1e-12);                 % out-of-fold GBM total
    pct_raw = 100 * contrib ./ Tg;             % % of GBM total
    row_sum = sum(pct_raw, 2);
    row_sum(row_sum <= 0) = 1e-12;

    pct_gbm  = pct_raw ./ (row_sum/100);       % rows sum to 100%
    frac_gbm = pct_gbm / 100;                  % 0–1 fractions

    FucoF   = frac_gbm(:,1);
    PerF    = frac_gbm(:,2);
    Hex19F  = frac_gbm(:,3);
    AlloF   = frac_gbm(:,4);
    But19F  = frac_gbm(:,5);
    ChlbF   = frac_gbm(:,6);
    ZeaF    = frac_gbm(:,7);
    DVChlF  = frac_gbm(:,8);

    % PSC
    PSC_micro = FucoF + PerF;
    PSC_nano  = Hex19F + But19F + AlloF;
    PSC_pico  = ChlbF + ZeaF + DVChlF;

    % PFTs
    PFT_diatoms       = FucoF;
    PFT_dino          = PerF;
    PFT_prym          = Hex19F + But19F;
    PFT_cryptophytes  = AlloF;
    PFT_prok_all      = ZeaF + DVChlF;
    PFT_picoeukGreen  = max(0, PSC_pico - PFT_prok_all);

    clamp01 = @(x) max(0, min(1, x));

    PSC_micro = clamp01(PSC_micro);
    PSC_nano  = clamp01(PSC_nano);
    PSC_pico  = clamp01(PSC_pico);

    PFT_diatoms      = clamp01(PFT_diatoms);
    PFT_dino         = clamp01(PFT_dino);
    PFT_prym         = clamp01(PFT_prym);
    PFT_cryptophytes = clamp01(PFT_cryptophytes);
    PFT_prok_all     = clamp01(PFT_prok_all);
    PFT_picoeukGreen = clamp01(PFT_picoeukGreen);

    % -------------------- Build bin-based empirical model --------------
    fracStruct = struct();
    fracStruct.PSC_micro       = PSC_micro;
    fracStruct.PSC_nano        = PSC_nano;
    fracStruct.PSC_pico        = PSC_pico;
    fracStruct.PFT_diatoms     = PFT_diatoms;
    fracStruct.PFT_dino        = PFT_dino;
    fracStruct.PFT_prym        = PFT_prym;
    fracStruct.PFT_cryptophytes= PFT_cryptophytes;
    fracStruct.PFT_prok_all    = PFT_prok_all;
    fracStruct.PFT_picoeukGreen= PFT_picoeukGreen;

    binModel = build_bin_model_from_training(Tchla, fracStruct, 100);

    % -------------------- Pack output struct ---------------------------
    model = struct();
    model.Ttrain    = T_use;
    model.Tchla     = Tchla;
    model.P         = P;
    model.nnlsCoeff = nnlsCoeff;
    model.TChla_cv  = TChla_cv;
    model.gbmFull   = gbmFull;

    model.PSC = struct('micro',PSC_micro, ...
                       'nano', PSC_nano, ...
                       'pico', PSC_pico);

    model.PFT = struct('diatoms',     PFT_diatoms, ...
                       'dino',        PFT_dino, ...
                       'prym',        PFT_prym, ...
                       'cryptophytes',PFT_cryptophytes, ...
                       'prok_all',    PFT_prok_all, ...
                       'picoeukGreen',PFT_picoeukGreen);

    model.binModel = binModel;
end
