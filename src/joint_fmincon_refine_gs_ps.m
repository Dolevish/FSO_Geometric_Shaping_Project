function out = joint_fmincon_refine_gs_ps(x0, p0, params, opts)
%JOINT_FMINCON_REFINE_GS_PS  Joint local refinement of constellation x and probabilities p.
%
% Solves the local optimization problem
%
%   maximize I(X;Y)
%
% over both x and p, under the IM/DD constraints:
%
%   x_i >= 0
%   x_{i+1} - x_i >= minGap
%   p_i >= p_min
%   sum_i p_i = 1
%   sum_i p_i x_i = P_avg
%
% The objective uses AMI_functions.AMI_noCSI_fast_grid. The returned solution
% is also validated with AMI_functions.AMI_noCSI_validate.
%
% Required files:
%   AMI_functions.m
%   calculate_Py_given_x.m
%   project_probabilities_power.m
%
% Optional files used for robust initial repair:
%   project_constellation_weighted_power.m
%
% Usage:
%   out = joint_fmincon_refine_gs_ps(x0, p0, params, opts);
%
% Main opts fields:
%   opts.minGap       default 0.005
%   opts.p_min        default 1e-8
%   opts.ghN_h        default 40
%   opts.xMaxUB       default max(8, 1.5*max(x0))
%   opts.nStarts      default 1
%   opts.xPerturbStd  default 0.03
%   opts.pPerturbStd  default 0.02
%   opts.maxIter      default 500
%   opts.maxFunEvals  default 20000
%   opts.display      default 'iter'
%   opts.seed         default 20260527

    if nargin < 4 || isempty(opts), opts = struct(); end

    x0 = real(x0(:));
    p0 = real(p0(:));
    M  = numel(x0);

    assert(numel(p0) == M, 'x0 and p0 must have the same length.');
    assert(isfield(params,'P_avg'), 'params.P_avg is required.');
    assert(isfield(params,'R'), 'params.R is required.');
    assert(isfield(params,'sigma_n_sq'), 'params.sigma_n_sq is required.');
    assert(isfield(params,'sigma_X_sq'), 'params.sigma_X_sq is required.');
    assert(isfield(params,'mu_t') && isfield(params,'sig_t'), 'params.mu_t and params.sig_t are required.');

    P_avg = params.P_avg;

    minGap      = get_opt(opts, 'minGap', 0.005);
    p_min       = get_opt(opts, 'p_min', 1e-8);
    ghN_h       = get_opt(opts, 'ghN_h', 40);
    nStarts     = get_opt(opts, 'nStarts', 1);
    xPerturbStd = get_opt(opts, 'xPerturbStd', 0.03);
    pPerturbStd = get_opt(opts, 'pPerturbStd', 0.02);
    maxIter     = get_opt(opts, 'maxIter', 500);
    maxFunEvals = get_opt(opts, 'maxFunEvals', 20000);
    displayMode = get_opt(opts, 'display', 'iter');
    seed        = get_opt(opts, 'seed', 20260527);
    useParallel = get_opt(opts, 'useParallelFiniteDiff', false);

    xMaxUB = get_opt(opts, 'xMaxUB', max(8, 1.5*max(x0)));
    xMaxUB = max(xMaxUB, max(x0) + 0.5);

    % Build a wide, constellation-independent grid for the objective.
    % Because fmincon may move x during the search, the grid is based on the
    % explicit upper bound rather than the initial maximum x only.
    y_grid = get_opt(opts, 'y_grid', []);
    if isempty(y_grid)
        y_grid = AMI_functions.build_noCSI_y_grid(params, xMaxUB);
    end

    % Feasibility check for minGap baseline: b=(0:M-1)*minGap must be able to
    % fit under weighted power for at least one legal p.
    if (M-1)*minGap > xMaxUB
        error('xMaxUB is too small for the requested minGap.');
    end

    % fmincon decision vector: u = [x; p]
    nvar = 2*M;

    % Linear inequality for minGap:
    % x_i - x_{i+1} <= -minGap
    A = zeros(M-1, nvar);
    b = -minGap * ones(M-1, 1);
    for i = 1:(M-1)
        A(i, i)   = 1;
        A(i, i+1) = -1;
    end

    % Linear equality for probability normalization: sum(p)=1
    Aeq = zeros(1, nvar);
    Aeq(M+1:end) = 1;
    beq = 1;

    lb = [zeros(M,1); p_min * ones(M,1)];
    ub = [xMaxUB * ones(M,1); ones(M,1)];

    obj = @(u) joint_objective(u, M, params, ghN_h, y_grid);
    nonlcon = @(u) joint_nonlcon(u, M, P_avg);

    fopts = optimoptions('fmincon', ...
        'Algorithm', 'sqp', ...
        'Display', displayMode, ...
        'MaxIterations', maxIter, ...
        'MaxFunctionEvaluations', maxFunEvals, ...
        'ConstraintTolerance', 1e-10, ...
        'OptimalityTolerance', 1e-7, ...
        'StepTolerance', 1e-10, ...
        'FiniteDifferenceType', 'central', ...
        'UseParallel', logical(useParallel));

    rng(seed, 'twister');

    % Build feasible starting points.
    U0 = cell(nStarts, 1);
    [x_start, p_start] = repair_xp(x0, p0, params, minGap, p_min, xMaxUB);
    U0{1} = [x_start; p_start];

    for s = 2:nStarts
        x_try = x_start + xPerturbStd * randn(M,1);
        p_try = p_start + pPerturbStd * randn(M,1);
        [x_try, p_try] = repair_xp(x_try, p_try, params, minGap, p_min, xMaxUB);
        U0{s} = [x_try; p_try];
    end

    results = repmat(struct( ...
        'startIndex', NaN, ...
        'u0', [], ...
        'u', [], ...
        'x', [], ...
        'p', [], ...
        'ami_fast', NaN, ...
        'ami_val', NaN, ...
        'exitflag', NaN, ...
        'output', [], ...
        'runtime', NaN, ...
        'constraint', []), nStarts, 1);

    fprintf('\n============================================================\n');
    fprintf('   Joint fmincon GS/PS refinement\n');
    fprintf('============================================================\n');
    fprintf('M=%d | P_avg=%.6f | minGap=%.6g | p_min=%.3e | xMaxUB=%.3f\n', ...
        M, P_avg, minGap, p_min, xMaxUB);
    fprintf('nStarts=%d | maxIter=%d | maxFunEvals=%d | ghN_h=%d | Ny=%d\n', ...
        nStarts, maxIter, maxFunEvals, ghN_h, numel(y_grid));

    for s = 1:nStarts
        fprintf('\n------------------------------------------------------------\n');
        fprintf(' Joint fmincon start %d/%d\n', s, nStarts);
        fprintf('------------------------------------------------------------\n');

        u0 = U0{s};
        x0s = u0(1:M);
        p0s = u0(M+1:end);
        ami0_fast = AMI_functions.AMI_noCSI_fast_grid(x0s(:).', p0s(:).', params, ghN_h, y_grid);
        ami0_val  = AMI_functions.AMI_noCSI_validate(x0s(:).', p0s(:).', params);

        fprintf('Initial x = ['); fprintf(' %.6f', x0s); fprintf(' ]\n');
        fprintf('Initial p = ['); fprintf(' %.6e', p0s); fprintf(' ]\n');
        fprintf('Initial constraints: sum(p)=%.12f | pTx=%.12f | minDiff=%.3e | min(p)=%.3e\n', ...
            sum(p0s), dot(p0s,x0s), min(diff(x0s)), min(p0s));
        fprintf('Initial AMI fast=%.8f | validated=%.8f\n', ami0_fast, ami0_val);

        t0 = tic;
        [u, fval, exitflag, output] = fmincon(obj, u0, A, b, Aeq, beq, lb, ub, nonlcon, fopts);
        runtime = toc(t0);

        x = u(1:M);
        p = u(M+1:end);

        % Numerical cleanup only; do not reorder after fmincon.
        x(abs(x) < 1e-12) = 0;
        p(p < p_min) = p_min;
        p = p / sum(p);
        % Re-project p only if numerical cleanup slightly moved pTx.
        if abs(dot(p,x) - P_avg) > 1e-8
            p = project_probabilities_power(p, x, P_avg, p_min);
        end
        u = [x; p];

        ami_fast = AMI_functions.AMI_noCSI_fast_grid(x(:).', p(:).', params, ghN_h, y_grid);
        ami_val  = AMI_functions.AMI_noCSI_validate(x(:).', p(:).', params);

        cinfo = constraint_info(x, p, P_avg, minGap);

        fprintf('\nResult start %d:\n', s);
        fprintf('    exitflag      = %d\n', exitflag);
        fprintf('    fval          = %.10f\n', fval);
        fprintf('    runtime       = %.2f sec\n', runtime);
        fprintf('    x             = ['); fprintf(' %.6f', x); fprintf(' ]\n');
        fprintf('    p             = ['); fprintf(' %.6e', p); fprintf(' ]\n');
        fprintf('    sum(p)        = %.12f\n', cinfo.sum_p);
        fprintf('    dot(p,x)      = %.12f\n', cinfo.power);
        fprintf('    min diff(x)   = %.12e\n', cinfo.min_diff);
        fprintf('    min(p)        = %.3e\n', cinfo.min_p);
        fprintf('    AMI fast      = %.8f bits/symbol\n', ami_fast);
        fprintf('    AMI validated = %.8f bits/symbol\n', ami_val);

        results(s).startIndex = s;
        results(s).u0         = u0;
        results(s).u          = u;
        results(s).x          = x;
        results(s).p          = p;
        results(s).ami_fast   = ami_fast;
        results(s).ami_val    = ami_val;
        results(s).exitflag   = exitflag;
        results(s).output     = output;
        results(s).runtime    = runtime;
        results(s).constraint = cinfo;
    end

    [bestVal, bestIdx] = max([results.ami_val]);

    out = struct();
    out.bestIndex   = bestIdx;
    out.x_best      = results(bestIdx).x;
    out.p_best      = results(bestIdx).p;
    out.ami_fast    = results(bestIdx).ami_fast;
    out.ami_val     = bestVal;
    out.exitflag    = results(bestIdx).exitflag;
    out.output      = results(bestIdx).output;
    out.results     = results;
    out.params      = params;
    out.opts        = opts;
    out.y_grid      = y_grid;

    fprintf('\n============================================================\n');
    fprintf(' Best joint fmincon result selected by validated AMI\n');
    fprintf('============================================================\n');
    fprintf('best start     = %d\n', bestIdx);
    fprintf('x_best         = ['); fprintf(' %.6f', out.x_best); fprintf(' ]\n');
    fprintf('p_best         = ['); fprintf(' %.6e', out.p_best); fprintf(' ]\n');
    fprintf('AMI fast       = %.8f bits/symbol\n', out.ami_fast);
    fprintf('AMI validated  = %.8f bits/symbol\n', out.ami_val);
    fprintf('sum(p)         = %.12f\n', sum(out.p_best));
    fprintf('dot(p,x)       = %.12f\n', dot(out.p_best, out.x_best));
    fprintf('min diff(x)    = %.12e\n', min(diff(out.x_best)));
    fprintf('min(p)         = %.3e\n', min(out.p_best));
end

% =========================================================================
% Objective and constraints
% =========================================================================
function f = joint_objective(u, M, params, ghN_h, y_grid)
    x = u(1:M);
    p = u(M+1:end);
    p = max(p, realmin);
    p = p / sum(p);

    mi = AMI_functions.AMI_noCSI_fast_grid(x(:).', p(:).', params, ghN_h, y_grid);
    if ~isfinite(mi)
        f = 1e6;
    else
        f = -mi;
    end
end

function [c, ceq] = joint_nonlcon(u, M, P_avg)
    x = u(1:M);
    p = u(M+1:end);
    c = [];
    ceq = dot(p, x) - P_avg;
end

% =========================================================================
% Feasible repair and reporting helpers
% =========================================================================
function [x, p] = repair_xp(x_in, p_in, params, minGap, p_min, xMaxUB)
    x = real(x_in(:));
    p = real(p_in(:));
    M = numel(x);

    p(~isfinite(p)) = p_min;
    p = max(p, p_min);
    p = p / sum(p);

    % Create minimal cfg object for project_constellation_weighted_power.
    cfg = struct();
    cfg.P_avg = params.P_avg;
    cfg.SA = struct();
    cfg.SA.minGap = minGap;

    if exist('project_constellation_weighted_power', 'file') == 2
        x = project_constellation_weighted_power(x, p, cfg);
    else
        x = sort(max(x, 0), 'ascend');
        b = (0:M-1).' * minGap;
        r = max(x - b, 0);
        for i = 2:M
            r(i) = max(r(i), r(i-1));
        end
        target = params.P_avg - dot(p, b);
        if target <= 0
            x = b;
        elseif dot(p,r) > 0
            x = b + r * (target / dot(p,r));
        else
            template = (0:M-1).';
            x = b + template * (target / dot(p, template));
        end
    end

    x = min(max(x, 0), xMaxUB);

    % If clipping changed the weighted power, project p for this fixed x.
    p = project_probabilities_power(p, x, params.P_avg, p_min);

    % Tiny numerical cleanup, no sorting after p projection.
    x(abs(x) < 1e-12) = 0;
    p(p < p_min) = p_min;
    p = p / sum(p);
    if abs(dot(p,x) - params.P_avg) > 1e-9
        p = project_probabilities_power(p, x, params.P_avg, p_min);
    end
end

function cinfo = constraint_info(x, p, P_avg, minGap)
    cinfo = struct();
    cinfo.sum_p = sum(p);
    cinfo.power = dot(p,x);
    cinfo.power_error = dot(p,x) - P_avg;
    cinfo.min_diff = min(diff(x));
    cinfo.min_gap_error = min(diff(x)) - minGap;
    cinfo.min_p = min(p);
end

function v = get_opt(s, name, default)
    if isfield(s, name) && ~isempty(s.(name))
        v = s.(name);
    else
        v = default;
    end
end
