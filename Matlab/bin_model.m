function binModel = bin_model(Tchla, fracStruct, nbin)
% BUILD_BIN_MODEL_FROM_TRAINING
%   Construct bin-based empirical relationships between PSC/PFT fractions
%   and log10(Tchla) using weighted 4th-degree polynomial fits.
%
% Inputs:
%   Tchla      - in situ total chlorophyll-a [mg m^-3], vector
%   fracStruct - struct with fields:
%                   PSC_micro, PSC_nano, PSC_pico,
%                   PFT_diatoms, PFT_dino, PFT_prym,
%                   PFT_cryptophytes, PFT_prok_all, PFT_picoeukGreen
%                each N x 1, fractions in [0,1]
%   nbin       - number of log-spaced bins (default 100)
%
% Output:
%   binModel   - struct with fields:
%       .labels   - cell array of variable names
%       .coeffs   - [V x 5] polynomial coefficients (a..e) for fractions
%       .edges    - bin edges (Tchla)
%       .centers  - bin centers (Tchla)
%       .stats    - struct array with weighted R^2, r, p_approx, nEff, nBinsUsed
%       .Tmin, .Tmax - min and max Tchla used

    if nargin < 3 || isempty(nbin)
        nbin = 100;
    end

    labels = { ...
        'PSC_micro','PSC_nano','PSC_pico', ...
        'PFT_diatoms','PFT_dino','PFT_prym', ...
        'PFT_cryptophytes','PFT_prok_all','PFT_picoeukGreen'};

    V = numel(labels);

    % Collect fractions into matrix (rows: variables, cols: samples)
    vars = cell(V,1);
    vars{1} = fracStruct.PSC_micro;
    vars{2} = fracStruct.PSC_nano;
    vars{3} = fracStruct.PSC_pico;
    vars{4} = fracStruct.PFT_diatoms;
    vars{5} = fracStruct.PFT_dino;
    vars{6} = fracStruct.PFT_prym;
    vars{7} = fracStruct.PFT_cryptophytes;
    vars{8} = fracStruct.PFT_prok_all;
    vars{9} = fracStruct.PFT_picoeukGreen;

    % -------------------- Log-spaced bins ------------------------------
    Tmin = min(Tchla(Tchla > 0));
    Tmax = max(Tchla);
    edges_log = linspace(log10(Tmin), log10(Tmax), nbin+1);
    edges     = 10.^edges_log;
    centers   = sqrt(edges(1:end-1).*edges(2:end));

    % Binned means/stds and counts
    M = zeros(V, nbin);
    S = zeros(V, nbin);
    N = zeros(1, nbin);

    for b = 1:nbin
        inb = (Tchla >= edges(b)) & (Tchla < edges(b+1));
        N(b) = sum(inb);
        for v = 1:V
            y = vars{v};
            if N(b) > 0
                M(v,b) = mean(y(inb), 'omitnan');
                S(v,b) = std( y(inb), 'omitnan');
            else
                M(v,b) = NaN;
                S(v,b) = NaN;
            end
        end
    end

    % -------------------- Weighted 4th-degree fits ---------------------
    coeffs = NaN(V, 5);
    stats  = struct('label',labels, ...
                    'R2w',NaN, 'r',NaN, 'p_approx',NaN, ...
                    'nEff',NaN, 'nBinsUsed',NaN);

    x_all = log10(centers(:));
    w_all = N(:);

    for v = 1:V
        y_all = M(v,:).';

        mask = isfinite(x_all) & isfinite(y_all) & (w_all > 0);
        x = x_all(mask);
        y = y_all(mask);
        w = w_all(mask);

        if numel(x) < 5
            % not enough bins; leave NaNs
            continue;
        end

        X = [ones(size(x)), x, x.^2, x.^3, x.^4];
        W = diag(w);

        beta = (X' * W * X) \ (X' * W * y);
        coeffs(v,:) = beta(:).';

        % stats
        yhat = X * beta;
        mu_w = sum(w .* y) / sum(w);
        SSEw = sum(w .* (y - yhat).^2);
        SSTw = sum(w .* (y - mu_w).^2);
        R2w  = 1 - SSEw / max(SSTw, eps);

        my   = mu_w;
        mhat = sum(w .* yhat) / sum(w);
        covw = sum(w .* (y - my) .* (yhat - mhat)) / sum(w);
        vy   = sum(w .* (y - my).^2) / sum(w);
        vhat = sum(w .* (yhat - mhat).^2) / sum(w);
        rw   = covw / sqrt(max(vy*vhat, eps));

        neff = (sum(w)^2) / sum(w.^2);
        dof  = max(neff - 2, 1);
        tstat = rw * sqrt(dof / max(1 - rw^2, eps));
        p_approx = 2 * tcdf(-abs(tstat), dof);

        stats(v).label     = labels{v};
        stats(v).R2w       = R2w;
        stats(v).r         = rw;
        stats(v).p_approx  = p_approx;
        stats(v).nEff      = neff;
        stats(v).nBinsUsed = numel(x);
    end

    % -------------------- Pack output ----------------------------------
    binModel = struct();
    binModel.labels  = labels;
    binModel.coeffs  = coeffs;  % fraction = a + b x + c x^2 + d x^3 + e x^4, x=log10(Tchla)
    binModel.edges   = edges;
    binModel.centers = centers;
    binModel.Tmin    = Tmin;
    binModel.Tmax    = Tmax;
    binModel.N       = N;
    binModel.M       = M;
    binModel.S       = S;
    binModel.stats   = stats;
end
