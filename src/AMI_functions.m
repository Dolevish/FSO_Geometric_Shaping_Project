classdef AMI_functions
%AMI_FUNCTIONS  Shared AMI and constellation-projection utilities.

methods(Static)

function mi_bits = AMI_noCSI_fast_grid(x, px, params, ghN_h, y_grid)
%AMI_NOCSI_FAST_GRID  Fast no-CSI AMI via 1D GH mixture + grid quadrature.

if nargin < 4 || isempty(ghN_h),  ghN_h = 40; end
if nargin < 5 || isempty(y_grid)
    y_grid = AMI_functions.build_noCSI_y_grid(params, max(abs(x)) * 1.5);
end

x  = x(:).';
px = px(:).';
M  = numel(x);
Ny = numel(y_grid);

% --- Special case: No turbulence (pure AWGN) ---
if params.sig_t < 1e-10
    means = params.R * x(:);
    d2    = bsxfun(@minus, y_grid, means).^2;
    logPyx = -0.5*log(2*pi*params.sigma_n_sq) - d2/(2*params.sigma_n_sq);

    logPx_Pyx = bsxfun(@plus, log(px(:)), logPyx);
    logPy = AMI_functions.logsumexp(logPx_Pyx, 1);

    MI = 0;
    for j = 1:M
        Pyx_j     = exp(logPyx(j,:));
        log_ratio = (logPyx(j,:) - logPy) / log(2);
        integrand = Pyx_j .* log_ratio;
        integrand(~isfinite(integrand)) = 0;
        MI = MI + px(j) * trapz(y_grid, integrand);
    end
    mi_bits = max(0, MI);
    return;
end

% --- General case: 1D GH over h + trapz over y ---
[z_h, w_h] = AMI_functions.gaussHermite(ghN_h);
W_h     = w_h(:) / sqrt(pi);
h_vals  = exp(params.mu_t + sqrt(2) * params.sig_t * z_h(:));

inv2s2  = 1 / (2 * params.sigma_n_sq);
logNorm = -0.5 * log(2 * pi * params.sigma_n_sq);
logW_h  = log(W_h);

logPyx = zeros(M, Ny);
for j = 1:M
    means_j  = params.R * h_vals * x(j);
    d2       = bsxfun(@minus, y_grid, means_j).^2;
    logTerms = bsxfun(@plus, logW_h + logNorm, -inv2s2 * d2);
    logPyx(j,:) = AMI_functions.logsumexp(logTerms, 1);
end

logPx_Pyx = bsxfun(@plus, log(px(:)), logPyx);
logPy     = AMI_functions.logsumexp(logPx_Pyx, 1);

MI = 0;
for j = 1:M
    Pyx_j     = exp(logPyx(j,:));
    log_ratio = (logPyx(j,:) - logPy) / log(2);
    integrand = Pyx_j .* log_ratio;
    integrand(~isfinite(integrand)) = 0;
    MI = MI + px(j) * trapz(y_grid, integrand);
end

mi_bits = max(0, MI);
end


function mi_bits = AMI_noCSI_validate(x, px, params)
%AMI_NOCSI_VALIDATE  High-precision no-CSI AMI via adaptive integral + trapz.
%   Calls standalone calculate_Py_given_x.m for result validation.

x  = x(:).';
px = px(:).';
M  = numel(x);

y_grid = AMI_functions.build_val_y_grid(params, x);
Ny     = numel(y_grid);

Pyx = zeros(M, Ny);
for i = 1:M
    Pyx(i, :) = calculate_Py_given_x(y_grid, x(i), params);
end

mi_bits = AMI_functions.calculate_mi_quadrature(Pyx, px, y_grid);
end


% =========================================================================
%  CONSTELLATION PROJECTION
% =========================================================================

function x = project_constellation_1D(x_in, cfg)
%PROJECT_CONSTELLATION_1D  Project a 1-D constellation onto the feasible set.
%
% For the IM/DD + average-intensity case used by this project, the
% projection simultaneously enforces:
%   1) sorted constellation levels,
%   2) non-negativity,
%   3) optional x(1)=0 anchoring (cfg.SA.pinZero),
%   4) adjacent spacing >= cfg.SA.minGap,
%   5) mean(x) = cfg.P_avg.
%
% The implementation is idempotent (up to floating-point precision):
%       P(P(x)) = P(x)
% which is important because an objective evaluator must not silently change
% the constellation being scored.
%
% For the mean-power IM/DD case, write
%       x_i = (i-1)*d_min + q_i,
% where q_i is non-negative and non-decreasing.  The fixed spacing skeleton
% consumes mean power ((M-1)*d_min)/2; only q is rescaled to use the remaining
% power budget.  Therefore power normalization cannot destroy the min-gap.

    if ~isfield(cfg,'SA')
        error('AMI_functions:MissingSAConfig', 'cfg.SA is required.');
    end
    if ~isfield(cfg,'P_avg') || ~(isscalar(cfg.P_avg) && isfinite(cfg.P_avg) && cfg.P_avg > 0)
        error('AMI_functions:BadAveragePower', 'cfg.P_avg must be a finite positive scalar.');
    end

    sa = cfg.SA;
    x  = real(x_in(:));

    if isempty(x) || any(~isfinite(x))
        error('AMI_functions:BadConstellation', ...
            'Constellation must be non-empty and contain only finite values.');
    end

    % Backward-compatible defaults for existing scripts.
    if ~isfield(sa,'minGap') || isempty(sa.minGap), sa.minGap = 0; end
    if ~isfield(sa,'enforce_sort') || isempty(sa.enforce_sort), sa.enforce_sort = true; end
    if ~isfield(sa,'imdd_mode') || isempty(sa.imdd_mode), sa.imdd_mode = false; end
    if ~isfield(sa,'enforce_power') || isempty(sa.enforce_power), sa.enforce_power = true; end
    if ~isfield(sa,'powerConstraint') || isempty(sa.powerConstraint), sa.powerConstraint = "meansquare"; end
    if ~isfield(sa,'pinZero') || isempty(sa.pinZero), sa.pinZero = logical(sa.imdd_mode); end

    minGap = max(0, double(sa.minGap));
    mode   = string(sa.powerConstraint);
    M      = numel(x);

    if sa.enforce_sort
        x = sort(x, 'ascend');
    elseif minGap > 0
        error('AMI_functions:MinGapNeedsSorting', ...
            'minGap > 0 requires cfg.SA.enforce_sort = true.');
    end

    % ---------------------------------------------------------------------
    % Project the IM/DD mean-intensity formulation used by the IEEE work.
    % ---------------------------------------------------------------------
    if sa.imdd_mode && sa.enforce_power && mode == "mean"
        if sa.pinZero && M == 1
            error('AMI_functions:InfeasiblePinnedSingleton', ...
                'A one-point constellation cannot satisfy pinZero=true and P_avg>0 simultaneously.');
        end

        % Non-negativity / optional zero anchoring.
        if sa.pinZero
            x = x - min(x);          % x(1)=0 after sorting
        else
            x = max(x, 0);
        end

        % Fixed min-gap skeleton b_i=(i-1)d.
        b = (0:M-1).' * minGap;
        skeletonMean = mean(b);
        residualMean = cfg.P_avg - skeletonMean;

        tolPower = 100 * eps(max(1, cfg.P_avg));
        if residualMean < -tolPower
            error('AMI_functions:InfeasibleMinGap', ...
                ['Infeasible constraints: minGap=%.6g with M=%d requires at least ' ...
                 'mean power %.6g, but P_avg=%.6g.'], ...
                minGap, M, skeletonMean, cfg.P_avg);
        end
        residualMean = max(0, residualMean);

        % Remove the deterministic spacing skeleton.  Projection of q onto
        % the cone q_1<=q_2<=... with q_i>=0 is implemented by a cumulative
        % maximum.  This is deterministic and idempotent for feasible inputs.
        q = x - b;
        q = max(q, 0);
        q = cummax(q);

        if sa.pinZero
            % Preserve x(1)=0 => q(1)=0.
            q = q - q(1);
            q = max(q, 0);
            q = cummax(q);
        end

        if residualMean <= tolPower
            q(:) = 0;
        else
            mq = mean(q);
            if mq <= eps
                % Degenerate input (all points collapsed).  Build a simple
                % feasible residual instead of dividing by zero.
                q = zeros(M,1);
                if sa.pinZero
                    q(2:end) = residualMean * M / (M-1);
                else
                    q(:) = residualMean;
                end
            else
                q = q * (residualMean / mq);
            end
        end

        x = b + q;

    else
        % -----------------------------------------------------------------
        % Legacy/general path retained for non-IM/DD or mean-square cases.
        % The IEEE revision uses the branch above (IM/DD + mean intensity).
        % -----------------------------------------------------------------
        if sa.imdd_mode
            if sa.pinZero
                x = x - min(x);
            else
                x = max(x,0);
            end
        end

        if minGap > 0
            if sa.enforce_sort, x = sort(x,'ascend'); end
            x = AMI_functions.enforce_min_gap_sorted(x, minGap);
            if sa.imdd_mode && sa.pinZero
                x = x - min(x);
            end
        end

        if sa.enforce_power
            x = Enforce_Power_Constraint(x, cfg.P_avg, mode);
        end
    end

    if sa.enforce_sort
        x = sort(x, 'ascend');
    end

    % Fail fast if the projection did not produce the requested feasible set.
    AMI_functions.assert_constellation_feasible(x, cfg, 1e-10);
end


function assert_constellation_feasible(x, cfg, tol)
%ASSERT_CONSTELLATION_FEASIBLE  Validate constraints without modifying x.
%
% Throws a descriptive error if a constellation violates the constraints
% encoded in cfg.  The function is intentionally side-effect free.

    if nargin < 3 || isempty(tol), tol = 1e-10; end
    x = real(x(:));

    if ~isfield(cfg,'SA'), error('AMI_functions:MissingSAConfig','cfg.SA is required.'); end
    sa = cfg.SA;

    if ~isfield(sa,'minGap') || isempty(sa.minGap), sa.minGap = 0; end
    if ~isfield(sa,'enforce_sort') || isempty(sa.enforce_sort), sa.enforce_sort = true; end
    if ~isfield(sa,'imdd_mode') || isempty(sa.imdd_mode), sa.imdd_mode = false; end
    if ~isfield(sa,'enforce_power') || isempty(sa.enforce_power), sa.enforce_power = true; end
    if ~isfield(sa,'powerConstraint') || isempty(sa.powerConstraint), sa.powerConstraint = "meansquare"; end
    if ~isfield(sa,'pinZero') || isempty(sa.pinZero), sa.pinZero = logical(sa.imdd_mode); end

    if isempty(x) || any(~isfinite(x))
        error('AMI_functions:NonFiniteConstellation','Constellation contains non-finite values.');
    end

    if isfield(cfg,'M') && ~isempty(cfg.M) && numel(x) ~= cfg.M
        error('AMI_functions:WrongCardinality', ...
            'Expected M=%d points but received %d.', cfg.M, numel(x));
    end

    scale = max([1; abs(x); abs(cfg.P_avg)]);
    atol = tol * max(scale);

    if sa.enforce_sort && any(diff(x) < -atol)
        error('AMI_functions:NotSorted','Constellation is not sorted.');
    end

    if sa.imdd_mode && min(x) < -atol
        error('AMI_functions:NegativeIntensity', ...
            'IM/DD constellation contains a negative intensity %.6g.', min(x));
    end

    if sa.imdd_mode && sa.pinZero && abs(x(1)) > atol
        error('AMI_functions:ZeroNotPinned', ...
            'pinZero=true but the first level is %.6g.', x(1));
    end

    minGap = max(0,double(sa.minGap));
    if numel(x) > 1 && minGap > 0 && min(diff(x)) < minGap - atol
        error('AMI_functions:MinGapViolation', ...
            'Minimum adjacent gap is %.6g but %.6g is required.', min(diff(x)), minGap);
    end

    if sa.enforce_power
        mode = string(sa.powerConstraint);
        switch mode
            case "mean"
                p = mean(x);
            case "meansquare"
                p = mean(x.^2);
            otherwise
                error('AMI_functions:BadPowerMode','Unknown power constraint mode: %s', mode);
        end

        if abs(p - cfg.P_avg) > atol
            error('AMI_functions:PowerConstraintViolation', ...
                'Power constraint mismatch: actual %.16g, target %.16g.', p, cfg.P_avg);
        end
    end
end


function x = enforce_min_gap_sorted(x, minGap)
%ENFORCE_MIN_GAP_SORTED  Push x(k) >= x(k-1)+minGap (assumes sorted input).
    x = x(:);
    for k = 2:numel(x)
        x(k) = max(x(k), x(k-1) + minGap);
    end
end


% =========================================================================
%  MI CALCULATION METHODS
% =========================================================================

function MI = calculate_mi_quadrature(Py_given_x_matrix, symbol_probs, y_range)
%CALCULATE_MI_QUADRATURE  AMI via trapezoidal rule on a pre-computed grid.
    Pyx   = max(Py_given_x_matrix, realmin);
    px    = symbol_probs(:).';
    Py    = px * Pyx;
    logPy = log(max(Py, realmin));
    logPyx= log(Pyx);

    M     = size(Pyx,1);
    integrand_values = Pyx .* (logPyx - logPy) / log(2);

    terms = zeros(M,1);
    if numel(y_range) > 1
        for m = 1:M
            terms(m) = trapz(y_range, integrand_values(m,:));
        end
    end

    MI = max(0, sum(symbol_probs(:) .* terms));
end


% =========================================================================
%  GRID BUILDERS
% =========================================================================

function y_grid = build_noCSI_y_grid(params, x_max_bound)
%BUILD_NOCSI_Y_GRID  Non-uniform y-grid for IM/DD no-CSI SA evaluator.
%   Constellation-independent — build once per turbulence level, reuse in SA.

    if nargin < 2 || isempty(x_max_bound), x_max_bound = 5; end

    sig_n = sqrt(params.sigma_n_sq);

    if params.sig_t > 1e-6
        h_75   = exp(params.mu_t + norminv(0.75)   * params.sig_t);
        h_999  = exp(params.mu_t + norminv(0.999)  * params.sig_t);
        h_9999 = exp(params.mu_t + norminv(0.9999) * params.sig_t);
    else
        h_75 = 1; h_999 = 1; h_9999 = 1;
    end

    R = params.R;
    y_min   = -8 * sig_n;
    y_focus = R * h_75   * x_max_bound + 8*sig_n;
    y_tail  = R * h_999  * x_max_bound + 5*sig_n;
    y_max   = R * h_9999 * x_max_bound + 10*sig_n;

    y_focus = max(y_focus, y_min + 20*sig_n);
    y_tail  = max(y_tail,  y_focus + 10*sig_n);
    y_max   = max(y_max,   y_tail + 10*sig_n);

    dy_fine   = max(sig_n / 4, 1e-6);
    dy_medium = max(sig_n, dy_fine * 4);
    dy_coarse = max(sig_n * 4, dy_medium * 4);

    g_dense  = y_min : dy_fine   : y_focus;
    g_medium = (y_focus + dy_medium) : dy_medium : y_tail;
    g_coarse = (y_tail  + dy_coarse) : dy_coarse : y_max;

    y_grid = unique([g_dense, g_medium, g_coarse]);
end


function y_grid = build_val_y_grid(params, x)
%BUILD_VAL_Y_GRID  Dense non-uniform y-grid for validation (wider, finer).

    R   = params.R;
    sig = sqrt(params.sigma_n_sq);
    xM  = max(abs(x));

    if params.sig_t > 1e-6
        h_hi = exp(params.mu_t + norminv(1-1e-9) * params.sig_t);
        if ~isfinite(h_hi), h_hi = exp(params.mu_t + 8*params.sig_t); end
    else
        h_hi = 1;
    end

    y_abs = R * xM * h_hi + 10*sig;

    yc   = 12*sig;
    dy_c = max(sig/20, 1e-4);
    y_core = -yc : dy_c : yc;

    N_mid = 800;  N_tail = 800;
    nonzero   = abs(y_core(y_core ~= 0));
    start_mid = max(min(nonzero), dy_c);
    y_mid_max = min(0.1*y_abs, max(30*sig, 5*yc));

    if start_mid < y_mid_max
        y_mid_pos = exp(linspace(log(start_mid), log(y_mid_max), N_mid));
    else
        y_mid_pos = [];
    end

    start_tail = max(y_mid_max + dy_c, start_mid*10);
    if start_tail < y_abs
        y_tail_pos = exp(linspace(log(start_tail), log(y_abs), N_tail));
    else
        y_tail_pos = [];
    end

    y_pos  = unique([y_mid_pos, y_tail_pos]);
    y_grid = sort(unique([-fliplr(y_pos), y_core, y_pos]));
end


% =========================================================================
%  SHARED HELPERS
% =========================================================================

function [x, w] = gaussHermite(n)
%GAUSSHERMITE  Nodes/weights for int exp(-t^2) f(t) dt (Golub-Welsch).
    i = (1:n-1)'; a = sqrt(i/2);
    J = diag(a,1) + diag(a,-1);
    [V, D] = eig(J);
    [x, idx] = sort(diag(D));
    V = V(:,idx);
    w = sqrt(pi) * (V(1,:).^2)';
    x = x(:)';
end


function s = logsumexp(A, dim)
%LOGSUMEXP  Numerically stable log(sum(exp(A), dim)).
    if nargin < 2, dim = 1; end
    maxA = max(A, [], dim);
    s = maxA + log(sum(exp(bsxfun(@minus, A, maxA)), dim));
end


end % methods(Static)
end % classdef
