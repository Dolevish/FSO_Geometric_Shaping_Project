function out = simulated_annealing_joint_gs_ps(x0, p0, params, opts)
%SIMULATED_ANNEALING_JOINT_GS_PS
% Joint simulated annealing engine for GS+PS.
%
% This optimizer perturbs both the constellation locations x and the
% probabilities p, then repairs the proposal so that the constraints hold:
%
%   x_i >= 0
%   x_{i+1} - x_i >= minGap
%   p_i >= p_min
%   sum_i p_i = 1
%   sum_i p_i x_i = P_avg
%
% The objective is AMI_functions.AMI_noCSI_fast_grid(x,p,...).
% The output also includes standard validated AMI.
%
% Required files:
%   AMI_functions.m
%   calculate_Py_given_x.m
%   project_probabilities_power.m
%   project_constellation_weighted_power.m
%
% Usage:
%   out = simulated_annealing_joint_gs_ps(x0, p0, params, opts);

    if nargin < 4 || isempty(opts), opts = struct(); end

    x0 = real(x0(:));
    p0 = real(p0(:));
    M  = numel(x0);

    if numel(p0) ~= M
        error('x0 and p0 must have the same length.');
    end
    if ~isfield(params,'P_avg') || ~isfield(params,'sigma_n_sq') || ~isfield(params,'R')
        error('params must include P_avg, sigma_n_sq, and R.');
    end

    P_avg = params.P_avg;

    minGap     = get_opt(opts, 'minGap', 0.005);
    p_min      = get_opt(opts, 'p_min', 1e-8);
    xMaxUB     = get_opt(opts, 'xMaxUB', max(8, 1.5*max(x0)));
    ghN_h      = get_opt(opts, 'ghN_h', 40);
    y_grid     = get_opt(opts, 'y_grid', []);
    maxIter    = get_opt(opts, 'maxIter', 2500);
    T0         = get_opt(opts, 'T0', 0.025);
    Tend       = get_opt(opts, 'Tend', 1e-4);
    xStep0     = get_opt(opts, 'xStep0', 0.25);
    pStep0     = get_opt(opts, 'pStep0', 0.08);
    xStepEnd   = get_opt(opts, 'xStepEnd', 0.01);
    pStepEnd   = get_opt(opts, 'pStepEnd', 0.004);
    seed       = get_opt(opts, 'seed', 20260607);
    verbose    = get_opt(opts, 'verbose', false);
    displayEvery = get_opt(opts, 'displayEvery', 500);

    if isempty(y_grid)
        y_grid = AMI_functions.build_noCSI_y_grid(params, xMaxUB);
    end

    rng(seed, 'twister');

    [x, p] = repair_joint_state(x0, p0, params, minGap, p_min, xMaxUB);
    mi = safe_fast_mi(x, p, params, ghN_h, y_grid);

    x_best = x;
    p_best = p;
    mi_best = mi;

    accepted = 0;
    improved = 0;

    if verbose
        fprintf('Joint-SA START | initMI=%.8f | maxIter=%d\n', mi, maxIter);
    end

    t0 = tic;
    for it = 1:maxIter
        frac = (it-1) / max(1, maxIter-1);
        T = T0 * (Tend / T0)^frac;
        xStep = xStep0 * (xStepEnd / xStep0)^frac;
        pStep = pStep0 * (pStepEnd / pStep0)^frac;

        % Joint proposal. A mild multiplicative component helps explore the
        % high-intensity tail; additive noise helps local refinement.
        x_prop = x + xStep * randn(M,1) + 0.10*xStep*(xMaxUB+1) * (rand(M,1)-0.5);
        p_prop = p + pStep * randn(M,1);

        [x_prop, p_prop] = repair_joint_state(x_prop, p_prop, params, minGap, p_min, xMaxUB);
        mi_prop = safe_fast_mi(x_prop, p_prop, params, ghN_h, y_grid);

        accept = false;
        if isfinite(mi_prop)
            d = mi_prop - mi;
            if d >= 0 || rand() < exp(d / max(T, realmin))
                accept = true;
            end
        end

        if accept
            x = x_prop;
            p = p_prop;
            mi = mi_prop;
            accepted = accepted + 1;

            if mi > mi_best
                x_best = x;
                p_best = p;
                mi_best = mi;
                improved = improved + 1;
            end
        end

        if verbose && (mod(it, displayEvery) == 0 || it == maxIter)
            fprintf('Joint-SA iter %6d/%6d | cur=%.8f | best=%.8f | acc=%.3f\n', ...
                it, maxIter, mi, mi_best, accepted/it);
        end
    end

    runtime = toc(t0);

    ami_val = AMI_functions.AMI_noCSI_validate(x_best(:).', p_best(:).', params);

    out = struct();
    out.x_best = x_best;
    out.p_best = p_best;
    out.ami_fast = mi_best;
    out.ami_val = ami_val;
    out.accepted = accepted;
    out.improved = improved;
    out.acceptanceRate = accepted / max(1,maxIter);
    out.runtime = runtime;
    out.opts = opts;
    out.params = params;
    out.y_grid = y_grid;
end

% =========================================================================
% Helpers
% =========================================================================
function mi = safe_fast_mi(x, p, params, ghN_h, y_grid)
    try
        mi = AMI_functions.AMI_noCSI_fast_grid(x(:).', p(:).', params, ghN_h, y_grid);
        if ~isfinite(mi), mi = -Inf; end
    catch
        mi = -Inf;
    end
end

function [x, p] = repair_joint_state(x_in, p_in, params, minGap, p_min, xMaxUB)
    x = real(x_in(:));
    p = real(p_in(:));
    M = numel(x);

    p(~isfinite(p)) = p_min;
    p = max(p, p_min);
    p = p / sum(p);

    cfgLocal = struct();
    cfgLocal.P_avg = params.P_avg;
    cfgLocal.SA = struct();
    cfgLocal.SA.minGap = minGap;

    % First project x for the current p.
    x = project_constellation_weighted_power(x, p, cfgLocal);

    % If the x projection exceeded the explicit upper bound, compress the
    % dynamic range and repair p for the fixed x. This keeps SA inside the
    % same feasible box used by fmincon.
    if max(x) > xMaxUB
        b = (0:M-1).' * minGap;
        r = max(x - b, 0);
        for i = 2:M
            r(i) = max(r(i), r(i-1));
        end
        if max(b + r) > xMaxUB
            scaleUB = min(1, (xMaxUB - max(b)) / max(max(r), eps));
            r = scaleUB * r;
            x = b + r;
        end
    end

    x = min(max(x, 0), xMaxUB);

    % Repair p for fixed x to hit the weighted power exactly.
    p = project_probabilities_power(p, x, params.P_avg, p_min);

    % Final safety.
    p = max(p, p_min);
    p = p / sum(p);
    if abs(dot(p,x) - params.P_avg) > 1e-8
        p = project_probabilities_power(p, x, params.P_avg, p_min);
    end

    x(abs(x) < 1e-12) = 0;
end

function v = get_opt(s, name, default)
    if isfield(s, name) && ~isempty(s.(name))
        v = s.(name);
    else
        v = default;
    end
end
