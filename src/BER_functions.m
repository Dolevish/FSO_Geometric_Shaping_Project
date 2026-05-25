classdef BER_functions
%BER_FUNCTIONS  Shared utilities for BER calculation and ML decision thresholds.
%
%   This class provides static methods for:
%   - Computing ML-optimal decision thresholds for No-CSI FSO receivers
%   - Calculating BER using Gauss-Hermite quadrature over log-normal fading
%   - Computing P(y|x) for the turbulent FSO channel
%
%   The key insight for FSO channels with multiplicative fading (y = R·h·x + n):
%   - Higher-amplitude symbols have wider PDF spread due to h·x coupling
%   - ML-optimal thresholds are shifted toward lower symbols (not midpoints)
%   - Thresholds are found where P(y|x_i) = P(y|x_j)

methods(Static)

%% ========================================================================
%  MAIN BER CALCULATION
% =========================================================================

function ber = calculate_BER_noCSI_ML(x, params)
%CALCULATE_BER_NOCSI_ML  BER for No-CSI receiver with ML-optimal thresholds.
%
%   Computes BER by:
%   1. Finding ML-optimal decision thresholds (PDF intersection points)
%   2. Integrating SER over log-normal fading using Gauss-Hermite quadrature
%
%   Inputs:
%       x      - Constellation points (will be sorted)
%       params - Struct with: R, sigma_n_sq, sigma_X_sq
%
%   Output:
%       ber    - Bit error rate (Gray coding assumed)

    x = sort(x(:));
    M = numel(x);
    
    R = params.R;
    sigma_n = sqrt(params.sigma_n_sq);
    
    % Log-normal parameters for h
    sig_t_sq = log(1 + params.sigma_X_sq);
    sig_t = sqrt(sig_t_sq);
    mu_t = -0.5 * sig_t_sq;  % ensures E[h] = 1
    
    % Step 1: Compute ML-optimal decision thresholds
    thresholds = BER_functions.find_ML_thresholds_all(x, params);
    
    % Step 2: Compute BER using Gauss-Hermite integration over fading
    if sig_t < 1e-10
        % No turbulence (AWGN only, h = 1 deterministically)
        h = 1;
        ser = BER_functions.compute_ser_given_h(x, h, R, sigma_n, M, thresholds);
    else
        % Gauss-Hermite quadrature for fading
        nGH = 50;
        [z_h, w_h] = BER_functions.gaussHermite(nGH);
        w_h = w_h / sqrt(pi);
        
        % h values: h = exp(mu_t + sqrt(2)*sig_t*z)
        h_vals = exp(mu_t + sqrt(2) * sig_t * z_h);
        
        % Integrate SER over h
        ser = 0;
        for g = 1:nGH
            h = h_vals(g);
            w = w_h(g);
            ser_h = BER_functions.compute_ser_given_h(x, h, R, sigma_n, M, thresholds);
            ser = ser + w * ser_h;
        end
    end
    
    % SER to BER approximation (for M-PAM with Gray coding)
    ber = ser / log2(M);
    ber = max(min(ber, 0.5), 1e-15);
end


%% ========================================================================
%  ML THRESHOLD FINDING
% =========================================================================

function thresholds = find_ML_thresholds_all(x, params)
%FIND_ML_THRESHOLDS_ALL  Find all ML-optimal decision thresholds for constellation.
%
%   Returns vector of M+1 thresholds: [-inf, y*_1, y*_2, ..., y*_{M-1}, inf]
%   where y*_k is the ML threshold between x(k) and x(k+1).

    x = sort(x(:));
    M = numel(x);
    R = params.R;
    
    thresholds = zeros(M+1, 1);
    thresholds(1) = -inf;
    thresholds(end) = inf;
    
    for k = 1:M-1
        xi = x(k);
        xj = x(k+1);
        
        % Initial guess: Euclidean midpoint
        y_mid = R * (xi + xj) / 2;
        
        if params.sigma_X_sq < 1e-12
            % No turbulence (AWGN): midpoint is optimal
            thresholds(k+1) = y_mid;
        else
            % Find ML threshold: where P(y|x_i) = P(y|x_j)
            thresholds(k+1) = BER_functions.find_ML_threshold(xi, xj, y_mid, params);
        end
    end
end


function y_opt = find_ML_threshold(xi, xj, y_init, params)
%FIND_ML_THRESHOLD  Find the ML-optimal decision threshold between xi and xj.
%
%   The threshold y* satisfies: P(y*|xi) = P(y*|xj)
%   Equivalently: log P(y*|xi) - log P(y*|xj) = 0
%
%   Uses fzero with bracket expansion for robustness.

    % Define the log-likelihood difference function
    fun = @(y) BER_functions.log_Py_given_x(y, xi, params) - ...
               BER_functions.log_Py_given_x(y, xj, params);
    
    % Search for bracket around initial guess
    sigma_n = sqrt(params.sigma_n_sq);
    dy = max(5 * sigma_n, 0.1 * abs(xj - xi));
    
    a = y_init - 3*dy;
    b = y_init + 3*dy;
    fa = fun(a);
    fb = fun(b);
    
    % Expand bracket if needed
    max_iter = 10;
    iter = 0;
    while sign(fa) == sign(fb) && iter < max_iter
        a = a - dy;
        b = b + dy;
        fa = fun(a);
        fb = fun(b);
        iter = iter + 1;
    end
    
    % Find zero crossing
    if sign(fa) ~= sign(fb)
        try
            y_opt = fzero(fun, [a, b]);
        catch
            y_opt = y_init;  % fallback to midpoint
        end
    else
        % If no bracket found, search for minimum |f(y)|
        y_search = linspace(y_init - 10*dy, y_init + 10*dy, 1001);
        f_vals = arrayfun(fun, y_search);
        [~, idx] = min(abs(f_vals));
        y_opt = y_search(idx);
    end
end


function [y_b_refined, app_gap] = refine_boundary_continuous(y_init, xi, xj, ...
        constellation_points, symbol_probs, params)
%REFINE_BOUNDARY_CONTINUOUS  Refine boundary using calculate_Py_given_x directly.
%
%   This is an alternative method that uses the standalone calculate_Py_given_x
%   function for higher precision validation.
%
%   Also returns the APP gap at the solution for verification.

    fun = @(yy) log(max(calculate_Py_given_x(yy, xi, params), realmin)) ...
              - log(max(calculate_Py_given_x(yy, xj, params), realmin));

    dy = max(5*sqrt(params.sigma_n_sq), 1e-3);
    a = y_init - 3*dy;
    b = y_init + 3*dy;
    fa = fun(a);
    fb = fun(b);
    k = 0;
    
    while sign(fa) == sign(fb) && k < 8
        a = a - dy;
        b = b + dy;
        fa = fun(a);
        fb = fun(b);
        k = k + 1;
    end

    if sign(fa) ~= sign(fb)
        y_b_refined = fzero(fun, [a b]);
    else
        y_b_refined = y_init;
    end

    % Compute APP gap at refined boundary
    p_all = arrayfun(@(xx) max(calculate_Py_given_x(y_b_refined, xx, params), realmin), ...
                     constellation_points);
    post = (symbol_probs(:)' .* p_all) / sum(symbol_probs(:)' .* p_all);
    idx_i = find(constellation_points == xi, 1);
    idx_j = find(constellation_points == xj, 1);
    
    if ~isempty(idx_i) && ~isempty(idx_j)
        app_gap = abs(post(idx_i) - post(idx_j));
    else
        app_gap = NaN;
    end
end


function y_b = find_boundary_loglike(y, Li, Lj, y_mid)
%FIND_BOUNDARY_LOGLIKE  Find boundary from pre-computed P(y|x) on a grid.
%
%   Uses interpolation to find where log P(y|x_i) = log P(y|x_j).
%   This is useful when P(y|x) has already been computed on a dense grid.
%
%   Inputs:
%       y     - y-grid (vector)
%       Li    - P(y|x_i) evaluated on grid
%       Lj    - P(y|x_j) evaluated on grid
%       y_mid - Initial guess (Euclidean midpoint)

    y = y(:).';
    Li = Li(:).';
    Lj = Lj(:).';
    logLdiff = log(max(Li, realmin)) - log(max(Lj, realmin));
    f = @(yy) interp1(y, logLdiff, yy, 'pchip');

    dy = y(2) - y(1);
    r = 10*dy;
    rmax = (y(end) - y(1)) / 2;
    a = max(y(1), y_mid - r);
    b = min(y(end), y_mid + r);
    fa = f(a);
    fb = f(b);
    
    while sign(fa) == sign(fb) && r < rmax
        r = 2*r;
        a = max(y(1), y_mid - r);
        b = min(y(end), y_mid + r);
        fa = f(a);
        fb = f(b);
    end

    if sign(fa) ~= sign(fb)
        y_b = fzero(f, [a, b]);
    else
        yd = linspace(y(1), y(end), 20001);
        [~, k] = min(abs(f(yd)));
        y_b = yd(k);
    end
end


%% ========================================================================
%  P(y|x) COMPUTATION
% =========================================================================

function logP = log_Py_given_x(y, x_sym, params)
%LOG_PY_GIVEN_X  Compute log P(y|x) for the FSO turbulent channel.
%
%   P(y|x) = ∫ P(y|x,h) · P(h) dh
%   where P(y|x,h) = N(y; R·h·x, σ_n²) and h ~ LogNormal(μ_t, σ_t²)
%
%   Uses Gauss-Hermite quadrature for efficient computation.
%   This is faster than calculate_Py_given_x but slightly less accurate.

    R = params.R;
    sigma_n_sq = params.sigma_n_sq;
    
    sig_t_sq = log(1 + params.sigma_X_sq);
    sig_t = sqrt(sig_t_sq);
    mu_t = -0.5 * sig_t_sq;
    
    % Handle AWGN case
    if sig_t < 1e-10
        h = 1;
        mean_y = R * h * x_sym;
        logP = -0.5*log(2*pi*sigma_n_sq) - (y - mean_y)^2 / (2*sigma_n_sq);
        return;
    end
    
    % Gauss-Hermite quadrature
    nGH = 30;
    [z_h, w_h] = BER_functions.gaussHermite(nGH);
    w_h = w_h / sqrt(pi);
    
    h_vals = exp(mu_t + sqrt(2) * sig_t * z_h);
    
    % Compute in log domain for stability
    log_terms = zeros(nGH, 1);
    for g = 1:nGH
        h = h_vals(g);
        mean_y = R * h * x_sym;
        log_gauss = -0.5*log(2*pi*sigma_n_sq) - (y - mean_y)^2 / (2*sigma_n_sq);
        log_terms(g) = log(w_h(g)) + log_gauss;
    end
    
    % Log-sum-exp for numerical stability
    max_log = max(log_terms);
    logP = max_log + log(sum(exp(log_terms - max_log)));
end


%% ========================================================================
%  SER COMPUTATION
% =========================================================================

function ser = compute_ser_given_h(x, h, R, sigma_n, M, thresholds)
%COMPUTE_SER_GIVEN_H  SER for given h with pre-computed thresholds.
%
%   Thresholds are pre-computed ML-optimal values (do NOT scale with h).
%   Only the received signal mean scales with h: mean_y = R · h · x_i

    ser = 0;
    for i = 1:M
        % Received signal mean depends on actual h
        mean_y = R * h * x(i);
        
        % Lower threshold error (detected as symbol i-1)
        if i > 1
            p_low = normcdf(thresholds(i), mean_y, sigma_n);
        else
            p_low = 0;
        end
        
        % Upper threshold error (detected as symbol i+1)
        if i < M
            p_high = 1 - normcdf(thresholds(i+1), mean_y, sigma_n);
        else
            p_high = 0;
        end
        
        ser = ser + (1/M) * (p_low + p_high);
    end
end


function [avg_ser, symbol_error_probs] = compute_SER_quadgk(constellation_points, ...
        symbol_probs, boundaries, params)
%COMPUTE_SER_QUADGK  SER via high-precision adaptive integration.
%
%   Integrates P(y|x) over decision regions using quadgk.
%   This is slower but more accurate than GH-based methods.
%   Useful for validation.
%
%   Inputs:
%       constellation_points - Symbol values
%       symbol_probs         - Prior probabilities
%       boundaries           - Decision boundaries [−inf, y*_1, ..., y*_{M-1}, inf]
%       params               - Channel parameters

    M = numel(constellation_points);
    symbol_error_probs = nan(1, M);
    
    for i = 1:M
        a = boundaries(i);
        if ~isfinite(a), a = -inf; end
        b = boundaries(i+1);
        if ~isfinite(b), b = inf; end
        
        integrand = @(yy) calculate_Py_given_x(yy, constellation_points(i), params);
        
        try
            Pc = quadgk(integrand, a, b, 'AbsTol', 1e-12, 'RelTol', 1e-9);
            symbol_error_probs(i) = 1 - max(0, min(1, Pc));
        catch ME
            warning('quadgk failed for symbol %d: %s', i, ME.message);
        end
    end
    
    avg_ser = sum(symbol_probs(:)' .* symbol_error_probs);
end


%% ========================================================================
%  GAUSS-HERMITE QUADRATURE
% =========================================================================

function [x, w] = gaussHermite(n)
%GAUSSHERMITE  Gauss-Hermite quadrature nodes and weights.
%
%   Computes nodes x_i and weights w_i for the integral:
%       ∫_{-∞}^{∞} exp(-t²) f(t) dt ≈ Σ w_i f(x_i)
%
%   Uses the Golub-Welsch algorithm (eigenvalue decomposition of
%   the symmetric tridiagonal Jacobi matrix).

    i = (1:n-1)';
    a = sqrt(i/2);
    J = diag(a, 1) + diag(a, -1);
    [V, D] = eig(J);
    [x, idx] = sort(diag(D));
    V = V(:, idx);
    w = sqrt(pi) * (V(1,:).^2)';
    x = x(:);
    w = w(:);
end


end % methods(Static)
end % classdef
