function stats = validate_model(T_val, model)
% VALIDATE_MODEL
%   Independent validation using OCChla (OC-CCI) and the bin-based model
%   versus GBM+NNLS-derived PSC/PFT concentrations on the validation subset.
%
% Inputs:
%   T_val - validation table from prepare_so_matchups (30% per year)
%   model - struct from train_model (contains nnlsCoeff, gbmFull, binModel)
%
% Output:
%   stats - struct array with fields:
%       .name   - variable name (e.g., 'PSC_micro')
%       .n      - number of matched samples
%       .r      - Pearson correlation (log10-space)
%       .R2     - r^2
%       .RMSE   - RMSE in log10 space
%       .MAE    - MAE in log10 space
%
% Notes:
%   - No figures are generated.
%   - Concentrations compared in mg m^-3, but metrics are in log10 space.

    needed = {'Tchla','OCChla','lat','Fuco','Per','X19hex','Allo','X19but','Chl_b','Zea','DVChla'};
    if ~all(ismember(needed, T_val.Properties.VariableNames))
        error('T_val must contain columns: %s', strjoin(needed, ', '));
    end

    % -------------------- Extract variables ----------------------------
    Tchla  = T_val.Tchla;
    OCChla = T_val.OCChla;
    lat    = T_val.lat;

    P = [T_val.Fuco,   T_val.Per,   T_val.X19hex, ...
         T_val.Allo,   T_val.X19but, T_val.Chl_b, ...
         T_val.Zea,    T_val.DVChla];

    % -------------------- Apply same cleaning rules --------------------
    mask_lat = isfinite(lat) & lat >= -90 & lat <= -40;
    mask_T   = isfinite(Tchla) & Tchla > 0 & Tchla < 15;

    P_meas = P;
    P_meas(~isfinite(P_meas)) = NaN;
    measured = P_meas >= 1e-3;
    has4 = sum(measured,2) >= 4;

    keep = mask_lat & mask_T & has4 & isfinite(OCChla) & OCChla > 0;

    Tchla  = Tchla(keep);
    OCChla = OCChla(keep);
    P      = P(keep,:);

    P(~isfinite(P)) = 0;
    P(P < 0)        = 0;
    P(P <= 1e-3)    = 0;

    fprintf('Validation samples retained after filters: %d\n', numel(Tchla));

    % -------------------- GBM+NNLS PSC/PFT (as in training) ------------
    C_nnls  = model.nnlsCoeff(:);       % [8 x 1]
    contrib = P .* C_nnls.';           % [N x 8]
    DP_lin  = sum(contrib, 2);

    epsP  = 1e-6;
    Xfeat = [P, log(P + epsP), DP_lin, log(DP_lin + epsP)];

    TChla_gbm = predict(model.gbmFull, Xfeat);
    Tg = max(TChla_gbm, 1e-12);

    pct_raw = 100 * contrib ./ Tg;
    row_sum = sum(pct_raw, 2);
    row_sum(row_sum <= 0) = 1e-12;

    pct_gbm  = pct_raw ./ (row_sum/100);
    frac_gbm = pct_gbm / 100;

    FucoF   = frac_gbm(:,1);
    PerF    = frac_gbm(:,2);
    Hex19F  = frac_gbm(:,3);
    AlloF   = frac_gbm(:,4);
    But19F  = frac_gbm(:,5);
    ChlbF   = frac_gbm(:,6);
    ZeaF    = frac_gbm(:,7);
    DVChlF  = frac_gbm(:,8);

    PSC_micro = FucoF + PerF;
    PSC_nano  = Hex19F + But19F + AlloF;
    PSC_pico  = ChlbF + ZeaF + DVChlF;

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

    % GBM concentrations (mg m^-3)
    GBM.PSC_micro = PSC_micro .* TChla_gbm;
    GBM.PSC_nano  = PSC_nano  .* TChla_gbm;
    GBM.PSC_pico  = PSC_pico  .* TChla_gbm;

    GBM.Diatoms        = PFT_diatoms      .* TChla_gbm;
    GBM.Dinoflagellates= PFT_dino         .* TChla_gbm;
    GBM.Prymnesiophytes= PFT_prym         .* TChla_gbm;
    GBM.Cryptophytes   = PFT_cryptophytes .* TChla_gbm;
    GBM.ProkProchloro  = PFT_prok_all     .* TChla_gbm;
    GBM.PicoEukGreen   = PFT_picoeukGreen .* TChla_gbm;

    % -------------------- Apply bin model to OCChla --------------------
    bm = model.binModel;
    labels = bm.labels;
    coeffs = bm.coeffs;

    x = log10(OCChla);
    fracCop = struct();

    for v = 1:numel(labels)
        beta = coeffs(v,:);
        a = beta(1); b = beta(2); c = beta(3); d = beta(4); e = beta(5);

        y = a + b*x + c*x.^2 + d*x.^3 + e*x.^4;   % fraction
        y = max(0, min(1, y));                   % clamp to 0–1

        fracCop.(labels{v}) = y;
    end

    Cop.PSC_micro = fracCop.PSC_micro .* OCChla;
    Cop.PSC_nano  = fracCop.PSC_nano  .* OCChla;
    Cop.PSC_pico  = fracCop.PSC_pico  .* OCChla;

    Cop.Diatoms        = fracCop.PFT_diatoms      .* OCChla;
    Cop.Dinoflagellates= fracCop.PFT_dino         .* OCChla;
    Cop.Prymnesiophytes= fracCop.PFT_prym         .* OCChla;
    Cop.Cryptophytes   = fracCop.PFT_cryptophytes .* OCChla;
    Cop.ProkProchloro  = fracCop.PFT_prok_all     .* OCChla;
    Cop.PicoEukGreen   = fracCop.PFT_picoeukGreen .* OCChla;

    % -------------------- Skill metrics in log10 space -----------------
    names = {'PSC_micro','PSC_nano','PSC_pico', ...
             'Diatoms','Dinoflagellates','Prymnesiophytes', ...
             'Cryptophytes','ProkProchloro','PicoEukGreen'};

    stats = struct('name',names, ...
                   'n',NaN,'r',NaN,'R2',NaN,'RMSE',NaN,'MAE',NaN);

    for i = 1:numel(names)
        g = GBM.(names{i});
        c = Cop.(names{i});

        m = isfinite(g) & isfinite(c) & g > 0 & c > 0;
        if nnz(m) < 10
            stats(i).n = nnz(m);
            continue;
        end

        lg = log10(g(m));
        lc = log10(c(m));

        R = corrcoef(lg, lc);
        r = R(1,2);
        R2 = r^2;
        RMSE = sqrt(mean((lc - lg).^2));
        MAE  = mean(abs(lc - lg));

        stats(i).n    = numel(lg);
        stats(i).r    = r;
        stats(i).R2   = R2;
        stats(i).RMSE = RMSE;
        stats(i).MAE  = MAE;
    end
end
