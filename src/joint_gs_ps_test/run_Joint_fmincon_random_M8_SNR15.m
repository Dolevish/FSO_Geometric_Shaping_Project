function run_Joint_fmincon_random_M8_SNR15()
%RUN_JOINT_FMINCON_RANDOM_M8_SNR15
% Joint GS+PS fmincon refinement starting ONLY from random feasible points.
%
% This script intentionally does NOT initialize from the Alternating GS/PS
% solution. Every fmincon start is generated randomly and then projected to
% satisfy the constraints:
%
%   x_i >= 0
%   x_{i+1}-x_i >= minGap
%   p_i >= p_min
%   sum_i p_i = 1
%   sum_i p_i x_i = P_avg
%
% The final winner is selected by validated AMI, not by the fast objective.
%
% Required files:
%   AMI_functions.m
%   calculate_Py_given_x.m
%   define_constellation.m
%   project_probabilities_power.m
%   project_constellation_weighted_power.m
%   joint_fmincon_refine_gs_ps.m

    clear; clc; close all;

    fprintf('============================================================\n');
    fprintf('   Joint fmincon GS/PS from RANDOM feasible starts\n');
    fprintf('============================================================\n\n');

    %% ------------------------------------------------------------
    % 1. Case parameters
    % -------------------------------------------------------------
    M            = 8;
    P_avg        = 1;
    SNR_dB       = 15;
    sigma_X_sq   = 0.1;
    R            = 1;
    minGap       = 0.005;
    p_min        = 1e-8;
    ghN_h        = 40;
    xMaxUB       = 8;

    params = struct();
    params.M          = M;
    params.P_avg      = P_avg;
    params.SNR_dB     = SNR_dB;
    params.sigma_X_sq = sigma_X_sq;
    params.sigma_n_sq = P_avg^2 / (10^(SNR_dB/10));
    params.R          = R;

    sig_t_sq     = log(1 + sigma_X_sq);
    params.sig_t = sqrt(sig_t_sq);
    params.mu_t  = -0.5 * sig_t_sq;

    x_pam = define_constellation(M, P_avg, true, "mean");
    p_uni = ones(M,1)/M;
    ami_pam_fast = AMI_functions.AMI_noCSI_fast_grid( ...
        x_pam(:).', p_uni(:).', params, ghN_h, AMI_functions.build_noCSI_y_grid(params, xMaxUB));
    ami_pam_val  = AMI_functions.AMI_noCSI_validate(x_pam(:).', p_uni(:).', params);

    y_grid = AMI_functions.build_noCSI_y_grid(params, xMaxUB);

    fprintf('[1] Case parameters\n');
    fprintf('    M              = %d\n', M);
    fprintf('    SNR            = %.1f dB\n', SNR_dB);
    fprintf('    P_avg          = %.6f\n', P_avg);
    fprintf('    sigma_X_sq     = %.6f\n', sigma_X_sq);
    fprintf('    sigma_n_sq     = %.6e\n', params.sigma_n_sq);
    fprintf('    minGap         = %.6f\n', minGap);
    fprintf('    p_min          = %.2e\n', p_min);
    fprintf('    xMaxUB         = %.3f\n', xMaxUB);
    fprintf('    y_grid range   = [%.4f, %.4f], Ny=%d\n', y_grid(1), y_grid(end), numel(y_grid));

    fprintf('\n[2] Uniform PAM baseline\n');
    print_solution('Uniform PAM', x_pam(:), p_uni(:), P_avg, minGap, ami_pam_fast, ami_pam_val, ami_pam_val);

    %% ------------------------------------------------------------
    % 2. Random-start settings
    % -------------------------------------------------------------
    nRandomStarts = 10;       % Increase to 20+ for a stronger global search.
    seed          = 20260607;

    % fmincon settings for EACH random start. We call the existing joint
    % function with nStarts=1, so every top-level start is independent and random.
    opts = struct();
    opts.minGap       = minGap;
    opts.p_min        = p_min;
    opts.ghN_h        = ghN_h;
    opts.xMaxUB       = xMaxUB;
    opts.y_grid       = y_grid;
    opts.nStarts      = 1;
    opts.xPerturbStd  = 0;    % ignored because nStarts=1
    opts.pPerturbStd  = 0;    % ignored because nStarts=1
    opts.maxIter      = 700;
    opts.maxFunEvals  = 40000;
    opts.display      = 'iter';
    opts.useParallelFiniteDiff = false;

    rng(seed, 'twister');

    fprintf('\n[3] Random feasible starts\n');
    fprintf('    nRandomStarts = %d\n', nRandomStarts);
    fprintf('    seed          = %d\n', seed);
    fprintf('    NOTE: no start is initialized from Alternating/GS/PS solution.\n');

    allRuns = repmat(struct( ...
        'startIndex', NaN, ...
        'x0', [], 'p0', [], ...
        'ami0_fast', NaN, 'ami0_val', NaN, ...
        'out', [], ...
        'x_best', [], 'p_best', [], ...
        'ami_fast', NaN, 'ami_val', NaN, ...
        'gain_val', NaN), nRandomStarts, 1);

    tAll = tic;

    for s = 1:nRandomStarts
        fprintf('\n############################################################\n');
        fprintf(' RANDOM START %d/%d\n', s, nRandomStarts);
        fprintf('############################################################\n');

        [x0, p0] = make_random_feasible_start(M, P_avg, minGap, p_min, xMaxUB, s);

        ami0_fast = AMI_functions.AMI_noCSI_fast_grid(x0(:).', p0(:).', params, ghN_h, y_grid);
        ami0_val  = AMI_functions.AMI_noCSI_validate(x0(:).', p0(:).', params);

        fprintf('\nRandom feasible initial point %d\n', s);
        print_solution('Random initial solution', x0, p0, P_avg, minGap, ami0_fast, ami0_val, ami_pam_val);

        % Give each local run its own deterministic seed. The local function will
        % not create extra starts because opts.nStarts=1.
        opts.seed = seed + 1000*s;

        out_s = joint_fmincon_refine_gs_ps(x0, p0, params, opts);

        x_best = out_s.x_best(:);
        p_best = out_s.p_best(:);
        ami_fast = out_s.ami_fast;
        ami_val  = out_s.ami_val;

        allRuns(s).startIndex = s;
        allRuns(s).x0 = x0;
        allRuns(s).p0 = p0;
        allRuns(s).ami0_fast = ami0_fast;
        allRuns(s).ami0_val = ami0_val;
        allRuns(s).out = out_s;
        allRuns(s).x_best = x_best;
        allRuns(s).p_best = p_best;
        allRuns(s).ami_fast = ami_fast;
        allRuns(s).ami_val = ami_val;
        allRuns(s).gain_val = ami_val - ami_pam_val;

        fprintf('\nCompleted random start %d:\n', s);
        print_solution('Best local result from this random start', ...
            x_best, p_best, P_avg, minGap, ami_fast, ami_val, ami_pam_val);
    end

    totalTime = toc(tAll);

    %% ------------------------------------------------------------
    % 3. Select final by validated AMI
    % -------------------------------------------------------------
    vals = [allRuns.ami_val];
    [bestVal, bestIdx] = max(vals);
    best = allRuns(bestIdx);

    fprintf('\n============================================================\n');
    fprintf(' Final summary: random-start joint fmincon\n');
    fprintf('============================================================\n');
    fprintf('-------------------------------------------------------------------------------------\n');
    fprintf(' idx | AMI init val | AMI final fast | AMI final val | Gain val | x_max | min p\n');
    fprintf('-------------------------------------------------------------------------------------\n');
    for s = 1:nRandomStarts
        fprintf(' %3d | %12.8f | %14.8f | %13.8f | %8.6f | %5.3f | %.1e\n', ...
            s, allRuns(s).ami0_val, allRuns(s).ami_fast, allRuns(s).ami_val, ...
            allRuns(s).gain_val, max(allRuns(s).x_best), min(allRuns(s).p_best));
    end
    fprintf('-------------------------------------------------------------------------------------\n');

    fprintf('\nBest random-start joint result selected by validated AMI:\n');
    fprintf('    best start    = %d\n', bestIdx);
    fprintf('    AMI val       = %.8f bits/symbol\n', bestVal);
    fprintf('    Gain val      = %.8f bits/symbol\n', bestVal - ami_pam_val);
    fprintf('    total runtime = %.1f min\n', totalTime/60);
    print_solution('Best random-start joint solution', ...
        best.x_best, best.p_best, P_avg, minGap, best.ami_fast, best.ami_val, ami_pam_val);

    %% ------------------------------------------------------------
    % 4. Save and plot
    % -------------------------------------------------------------
    result = struct();
    result.params        = params;
    result.opts          = opts;
    result.seed          = seed;
    result.nRandomStarts = nRandomStarts;
    result.x_pam         = x_pam(:);
    result.p_uni         = p_uni(:);
    result.ami_pam_fast  = ami_pam_fast;
    result.ami_pam_val   = ami_pam_val;
    result.allRuns       = allRuns;
    result.bestIdx       = bestIdx;
    result.x_best        = best.x_best;
    result.p_best        = best.p_best;
    result.ami_best_fast = best.ami_fast;
    result.ami_best_val  = best.ami_val;
    result.gain_best_val = best.ami_val - ami_pam_val;
    result.totalTime     = totalTime;

    save('result_Joint_fmincon_random_M8_SNR15.mat', 'result');
    fprintf('\n[Save] Saved result_Joint_fmincon_random_M8_SNR15.mat\n');

    fig1 = figure('Color','w','Name','Random-start joint fmincon summary');
    plot(1:nRandomStarts, [allRuns.ami0_val], 'o-', 'LineWidth', 1.4); hold on;
    plot(1:nRandomStarts, [allRuns.ami_val], 'x-', 'LineWidth', 1.4);
    yline(ami_pam_val, 'k:', 'LineWidth', 1.0, 'HandleVisibility','off');
    grid on; box on;
    xlabel('Random start index');
    ylabel('Validated AMI [bits/symbol]');
    title('Joint fmincon from random feasible starts');
    legend('Initial random feasible','After joint fmincon','Location','best');
    exportgraphics(fig1, 'Joint_fmincon_random_AMI_by_start_M8_SNR15.png', 'Resolution', 300);
    fprintf('[Plot] Saved Joint_fmincon_random_AMI_by_start_M8_SNR15.png\n');

    fig2 = figure('Color','w','Name','Best random-start joint x/p');
    tiledlayout(2,1);
    nexttile;
    stem(1:M, best.x_best, 'filled', 'LineWidth', 1.3);
    grid on; box on;
    xlabel('Symbol index'); ylabel('x_i');
    title(sprintf('Best random-start joint constellation, AMI_{val}=%.6f', best.ami_val));

    nexttile;
    stem(1:M, best.p_best, 'filled', 'LineWidth', 1.3);
    grid on; box on;
    xlabel('Symbol index'); ylabel('p_i');
    title('Best random-start joint probabilities');
    exportgraphics(fig2, 'Joint_fmincon_random_best_x_p_M8_SNR15.png', 'Resolution', 300);
    fprintf('[Plot] Saved Joint_fmincon_random_best_x_p_M8_SNR15.png\n');

    fprintf('\nDone.\n');
end

% =========================================================================
% Local helpers
% =========================================================================
function [x, p] = make_random_feasible_start(M, P_avg, minGap, p_min, xMaxUB, startIndex)
    maxAttempts = 500;

    for attempt = 1:maxAttempts
        % Random probability vector. The mixture of profiles gives both dense
        % and sparse-ish probability shapes while remaining fully random.
        alphaChoices = [0.35, 0.7, 1.0, 2.0];
        alpha = alphaChoices(1 + mod(startIndex + attempt - 1, numel(alphaChoices)));
        v = rand(M,1).^(1/alpha);
        p = p_min + (1 - M*p_min) * v / sum(v);

        % Random sorted nonnegative x. Force x(1)=0 and make sure some point is
        % above P_avg so that weighted-power projection is well behaved.
        span = P_avg + (xMaxUB - P_avg) * rand();
        raw = sort(span * rand(M,1).^1.5, 'ascend');
        raw(1) = 0;
        raw(end) = max(raw(end), P_avg + minGap*(M-1));
        raw = sort(raw, 'ascend');

        % Project x for fixed p under weighted average power and minGap.
        cfg = struct();
        cfg.P_avg = P_avg;
        cfg.SA = struct();
        cfg.SA.minGap = minGap;
        cfg.SA.enforce_sort = true;
        cfg.SA.imdd_mode = true;
        cfg.SA.enforce_power = true;
        cfg.SA.powerConstraint = "mean";

        try
            x = project_constellation_weighted_power(raw, p, cfg);
            x = x(:);
            p = p(:);

            % Numerical repair of p only if needed.
            p = project_probabilities_power(p, x, P_avg, p_min);

            ok = all(isfinite(x)) && all(isfinite(p)) && ...
                 all(x >= -1e-10) && all(p >= p_min - 1e-12) && ...
                 abs(sum(p)-1) < 1e-9 && abs(dot(p,x)-P_avg) < 1e-8 && ...
                 min(diff(x)) >= minGap - 1e-10 && max(x) <= xMaxUB + 1e-8;

            if ok
                x(abs(x) < 1e-12) = 0;
                return;
            end
        catch
            % Try again.
        end
    end

    error('Could not generate a random feasible start after %d attempts.', maxAttempts);
end

function print_solution(name, x, p, P_avg, minGap, ami_fast, ami_val, ami_pam_val)
    x = x(:); p = p(:);
    fprintf('\n    %s\n', name);
    fprintf('    x = ['); fprintf(' %.6f', x); fprintf(' ]\n');
    fprintf('    p = ['); fprintf(' %.6e', p); fprintf(' ]\n');
    fprintf('    sum(p)      = %.12f\n', sum(p));
    fprintf('    dot(p,x)    = %.12f\n', dot(p,x));
    fprintf('    target power= %.12f\n', P_avg);
    fprintf('    min(p)      = %.3e\n', min(p));
    fprintf('    min diff(x) = %.3e\n', min(diff(x)));
    fprintf('    minGap      = %.3e\n', minGap);
    fprintf('    x_max       = %.6f\n', max(x));
    fprintf('    AMI fast    = %.8f bits/symbol\n', ami_fast);
    fprintf('    AMI val     = %.8f bits/symbol\n', ami_val);
    fprintf('    gain val    = %.8f bits/symbol\n', ami_val - ami_pam_val);
end
