function run_Joint_SA_fmincon_sweep_parallel()
%RUN_JOINT_SA_FMINCON_SWEEP_PARALLEL
% Comprehensive parallel sweep for hybrid Joint-SA + Joint-fmincon GS/PS optimization.
%
% For every pair:
%   sigma_X_sq in [0, 0.1, 0.2, 0.3]
%   SNR_dB     in 0:5:30
%
% the script starts from RANDOM FEASIBLE constellations/probabilities and runs
% joint fmincon refinement over both x and p:
%
%   maximize I(X;Y)
%
% subject to:
%   x_i >= 0
%   x_{i+1} - x_i >= minGap
%   p_i >= p_min
%   sum_i p_i = 1
%   sum_i p_i x_i = P_avg
%
% For each random start:
%   random feasible start -> joint SA over x,p -> joint fmincon over x,p -> validation
%
% The winner for each case is selected by standard validated AMI, not by fast
% AMI. Optional strict validation is then applied to the final winner.
%
% Required files in the same MATLAB folder:
%   AMI_functions.m
%   calculate_Py_given_x.m
%   define_constellation.m
%   project_probabilities_power.m
%   project_constellation_weighted_power.m
%   joint_fmincon_refine_gs_ps.m
%   simulated_annealing_joint_gs_ps.m
%
% Output folder:
%   results_joint_SA_fmincon_sweep_parallel/
%
% Main outputs:
%   result_Joint_SA_fmincon_sweep_parallel.mat
%   result_Joint_SA_fmincon_sweep_parallel.csv
%   Joint_SA_fmincon_sweep_AMI_combined.png
%   Joint_SA_fmincon_sweep_gain_combined.png
%   Joint_SA_fmincon_sweep_AMI_by_turbulence_sigma_*.png
%   Joint_SA_fmincon_sweep_gain_by_turbulence_sigma_*.png

    clear; clc; close all;
    addpath(pwd);

    fprintf('============================================================\n');
    fprintf('   Parallel Hybrid Joint-SA + Joint-fmincon GS/PS sweep\n');
    fprintf('============================================================\n\n');

    %% ============================================================
    %  1. Ordered case grid
    % =============================================================
    turbulenceList = [ 0.1, 0.2, 0.3];
    snrList        = 0:5:30;

    [Tgrid, Sgrid] = ndgrid(turbulenceList, snrList);
    caseSigma = Tgrid(:);
    caseSNR   = Sgrid(:);
    nCases    = numel(caseSigma);

    %% ============================================================
    %  2. Global simulation parameters
    % =============================================================
    cfg = struct();

    % Channel / constellation
    cfg.M            = 8;
    cfg.P_avg        = 1;
    cfg.R            = 1;
    cfg.minGap       = 0.005;
    cfg.p_min        = 1e-8;
    cfg.xMaxUB       = 8;

    % Fast objective used inside fmincon
    cfg.ghN_h        = 40;

    % Random feasible starts per case
    cfg.nRandomStarts = 10;       % Increase to 20+ for stronger global search.
    cfg.baseSeed      = 20260607;

    % Joint Simulated Annealing settings per random start.
    % This is the hybrid warm-start stage before fmincon.
    cfg.useJointSA    = true;
    cfg.saMaxIter     = 3000;
    cfg.saT0          = 0.030;
    cfg.saTend        = 1e-4;
    cfg.saXStep0      = 0.30;
    cfg.saPStep0      = 0.10;
    cfg.saXStepEnd    = 0.015;
    cfg.saPStepEnd    = 0.004;

    % fmincon settings per random start
    cfg.maxIter       = 700;
    cfg.maxFunEvals   = 40000;
    cfg.display       = 'off';    % Keep output readable in parallel sweep.

    % Parallelization
    cfg.useParallelCases = true;  % Parallel over (sigma_X_sq, SNR) cases.
    cfg.poolProfile      = 'Processes';
    cfg.numWorkers       = [];    % [] lets MATLAB choose default.

    % Important: keep this false because the OUTER case loop is already parallel.
    cfg.useParallelFiniteDiff = false;

    % Validation
    cfg.doStrictValidation = true;
    cfg.strictConfigs = make_strict_configs();

    % Output
    outDir = fullfile(pwd, 'results_joint_SA_fmincon_sweep_parallel');
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    fprintf('[1] Case grid\n');
    fprintf('    turbulenceList  = ['); fprintf(' %.3g', turbulenceList); fprintf(' ]\n');
    fprintf('    snrList         = ['); fprintf(' %.1f', snrList); fprintf(' ] dB\n');
    fprintf('    total cases     = %d\n', nCases);

    fprintf('\n[2] Main settings\n');
    fprintf('    M               = %d\n', cfg.M);
    fprintf('    P_avg           = %.6f\n', cfg.P_avg);
    fprintf('    minGap          = %.6f\n', cfg.minGap);
    fprintf('    p_min           = %.2e\n', cfg.p_min);
    fprintf('    xMaxUB          = %.3f\n', cfg.xMaxUB);
    fprintf('    nRandomStarts   = %d\n', cfg.nRandomStarts);
    fprintf('    joint SA enabled= %d\n', cfg.useJointSA);
    fprintf('    joint SA maxIter= %d\n', cfg.saMaxIter);
    fprintf('    fmincon maxIter = %d\n', cfg.maxIter);
    fprintf('    ghN_h objective = %d\n', cfg.ghN_h);
    fprintf('    strict validation = %d\n', cfg.doStrictValidation);

    %% ============================================================
    %  3. Start parallel pool
    % =============================================================
    if cfg.useParallelCases
        ppool = gcp('nocreate');
        if isempty(ppool)
            fprintf('\n[3] Starting parallel pool...\n');
            if isempty(cfg.numWorkers)
                parpool(cfg.poolProfile);
            else
                parpool(cfg.poolProfile, cfg.numWorkers);
            end
        else
            fprintf('\n[3] Parallel pool already active: workers=%d\n', ppool.NumWorkers);
        end
    else
        fprintf('\n[3] Running serially because cfg.useParallelCases=false\n');
    end

    %% ============================================================
    %  4. Run all cases
    % =============================================================
    results = repmat(empty_case_result(cfg.M), nCases, 1);

    fprintf('\n[4] Running sweep...\n');
    tAll = tic;

    if cfg.useParallelCases
        parfor caseIdx = 1:nCases
            results(caseIdx) = run_one_case(caseIdx, caseSigma(caseIdx), caseSNR(caseIdx), cfg, outDir);
        end
    else
        for caseIdx = 1:nCases
            results(caseIdx) = run_one_case(caseIdx, caseSigma(caseIdx), caseSNR(caseIdx), cfg, outDir);
        end
    end

    totalRuntime = toc(tAll);

    %% ============================================================
    %  5. Build summary table and save
    % =============================================================
    summaryTable = build_summary_table(results, cfg.M);

    matFile = fullfile(outDir, 'result_Joint_SA_fmincon_sweep_parallel.mat');
    csvFile = fullfile(outDir, 'result_Joint_SA_fmincon_sweep_parallel.csv');
    save(matFile, 'results', 'summaryTable', 'cfg', 'turbulenceList', 'snrList', 'totalRuntime');
    writetable(summaryTable, csvFile);

    fprintf('\n============================================================\n');
    fprintf(' Final hybrid sweep summary\n');
    fprintf('============================================================\n');
    fprintf('    total runtime = %.2f min\n', totalRuntime/60);
    fprintf('    saved MAT     = %s\n', matFile);
    fprintf('    saved CSV     = %s\n', csvFile);

    disp(summaryTable(:, {'caseIndex','sigma_X_sq','SNR_dB','AMI_PAM_val','AMI_best_SA_val','AMI_joint_val','AMI_joint_strict_min','Gain_strict_min','bestStart','x_max'}));

    %% ============================================================
    %  6. Plot results
    % =============================================================
    make_all_plots(summaryTable, turbulenceList, snrList, outDir);

    fprintf('\nDone.\n');
end

% =========================================================================
% Per-case optimization
% =========================================================================
function res = run_one_case(caseIdx, sigma_X_sq, SNR_dB, cfg, outDir)
    M = cfg.M;
    P_avg = cfg.P_avg;

    res = empty_case_result(M);
    res.caseIndex = caseIdx;
    res.sigma_X_sq = sigma_X_sq;
    res.SNR_dB = SNR_dB;

    fprintf('[Case %02d] START | sigma_X_sq=%.3f | SNR=%.1f dB\n', caseIdx, sigma_X_sq, SNR_dB);
    tCase = tic;

    params = make_params(cfg, sigma_X_sq, SNR_dB);

    % Baseline Uniform PAM.
    x_pam = define_constellation(M, P_avg, true, "mean");
    p_uni = ones(M,1)/M;
    ami_pam_val = AMI_functions.AMI_noCSI_validate(x_pam(:).', p_uni(:).', params);

    y_grid = AMI_functions.build_noCSI_y_grid(params, cfg.xMaxUB);
    ami_pam_fast = AMI_functions.AMI_noCSI_fast_grid(x_pam(:).', p_uni(:).', params, cfg.ghN_h, y_grid);

    % fmincon options passed to the existing joint refinement function.
    opts = struct();
    opts.minGap       = cfg.minGap;
    opts.p_min        = cfg.p_min;
    opts.ghN_h        = cfg.ghN_h;
    opts.xMaxUB       = cfg.xMaxUB;
    opts.y_grid       = y_grid;
    opts.nStarts      = 1;
    opts.xPerturbStd  = 0;
    opts.pPerturbStd  = 0;
    opts.maxIter      = cfg.maxIter;
    opts.maxFunEvals  = cfg.maxFunEvals;
    opts.display      = cfg.display;
    opts.useParallelFiniteDiff = cfg.useParallelFiniteDiff;

    runList = repmat(empty_start_result(M), cfg.nRandomStarts, 1);

    for s = 1:cfg.nRandomStarts
        seedStart = cfg.baseSeed + 100000*caseIdx + 1000*s;
        [x0, p0] = make_random_feasible_start(M, P_avg, cfg.minGap, cfg.p_min, cfg.xMaxUB, seedStart);

        ami0_val = AMI_functions.AMI_noCSI_validate(x0(:).', p0(:).', params);
        ami0_fast = AMI_functions.AMI_noCSI_fast_grid(x0(:).', p0(:).', params, cfg.ghN_h, y_grid);

        % Hybrid warm-start: Joint SA over both x and p.
        x_sa = x0;
        p_sa = p0;
        ami_sa_fast = ami0_fast;
        ami_sa_val  = ami0_val;
        saRuntime   = 0;
        saAcceptRate = NaN;

        try
            if cfg.useJointSA
                saOpts = struct();
                saOpts.minGap       = cfg.minGap;
                saOpts.p_min        = cfg.p_min;
                saOpts.xMaxUB       = cfg.xMaxUB;
                saOpts.ghN_h        = cfg.ghN_h;
                saOpts.y_grid       = y_grid;
                saOpts.maxIter      = cfg.saMaxIter;
                saOpts.T0           = cfg.saT0;
                saOpts.Tend         = cfg.saTend;
                saOpts.xStep0       = cfg.saXStep0;
                saOpts.pStep0       = cfg.saPStep0;
                saOpts.xStepEnd     = cfg.saXStepEnd;
                saOpts.pStepEnd     = cfg.saPStepEnd;
                saOpts.seed         = seedStart + 11;
                saOpts.verbose      = false;

                evalc('out_sa = simulated_annealing_joint_gs_ps(x0, p0, params, saOpts);');
                x_sa = out_sa.x_best(:);
                p_sa = out_sa.p_best(:);
                ami_sa_fast = out_sa.ami_fast;
                ami_sa_val  = out_sa.ami_val;
                saRuntime   = out_sa.runtime;
                saAcceptRate = out_sa.acceptanceRate;
            end

            % fmincon starts from the SA-refined state, not from the raw random state.
            opts.seed = seedStart + 17;
            evalc('out_s = joint_fmincon_refine_gs_ps(x_sa, p_sa, params, opts);');
            x_best = out_s.x_best(:);
            p_best = out_s.p_best(:);
            ami_fast = out_s.ami_fast;
            ami_val = out_s.ami_val;
            exitflag = out_s.exitflag;
            ok = true;
            errMsg = '';
        catch ME
            x_best = nan(M,1);
            p_best = nan(M,1);
            ami_fast = -Inf;
            ami_val = -Inf;
            exitflag = NaN;
            ok = false;
            errMsg = ME.message;
        end

        runList(s).startIndex = s;
        runList(s).seed = seedStart;
        runList(s).x0 = x0;
        runList(s).p0 = p0;
        runList(s).ami0_fast = ami0_fast;
        runList(s).ami0_val = ami0_val;
        runList(s).x_sa = x_sa;
        runList(s).p_sa = p_sa;
        runList(s).ami_sa_fast = ami_sa_fast;
        runList(s).ami_sa_val = ami_sa_val;
        runList(s).saRuntime = saRuntime;
        runList(s).saAcceptRate = saAcceptRate;
        runList(s).x_best = x_best;
        runList(s).p_best = p_best;
        runList(s).ami_fast = ami_fast;
        runList(s).ami_val = ami_val;
        runList(s).exitflag = exitflag;
        runList(s).ok = ok;
        runList(s).errorMessage = errMsg;
    end

    vals = [runList.ami_val];
    [bestVal, bestIdx] = max(vals);
    bestRun = runList(bestIdx);

    xj = bestRun.x_best(:);
    pj = bestRun.p_best(:);

    % Strict validation of the winning solution.
    strictSummary = empty_strict_summary();
    if cfg.doStrictValidation && all(isfinite(xj)) && all(isfinite(pj))
        try
            strictSummary = strict_validate_solution(xj, pj, params, x_pam(:), p_uni(:), cfg.strictConfigs);
        catch ME
            strictSummary.errorMessage = ME.message;
        end
    end

    res.params = params;
    res.x_pam = x_pam(:);
    res.p_uni = p_uni(:);
    res.ami_pam_fast = ami_pam_fast;
    res.ami_pam_val = ami_pam_val;
    res.allStarts = runList;
    res.bestStart = bestIdx;
    res.x_best = xj;
    res.p_best = pj;
    res.ami_fast = bestRun.ami_fast;
    res.ami_val = bestVal;
    res.gain_val = bestVal - ami_pam_val;
    res.strict = strictSummary;
    res.runtime = toc(tCase);

    % Save each case separately so a long run can be inspected even before the
    % full sweep completes.
    safeSigma = strrep(sprintf('%.3f', sigma_X_sq), '.', 'p');
    caseFile = fullfile(outDir, sprintf('case_%02d_sigma_%s_SNR_%02d.mat', caseIdx, safeSigma, round(SNR_dB)));
    save_case_result(caseFile, res);

    fprintf('[Case %02d] DONE  | sigma=%.3f | SNR=%.1f | AMI_val=%.8f | strict_min=%.8f | bestStart=%d | %.1fs\n', ...
        caseIdx, sigma_X_sq, SNR_dB, res.ami_val, res.strict.ami_min, bestIdx, res.runtime);
end


function save_case_result(caseFile, res)
    result = res; %#ok<NASGU>
    save(caseFile, 'result');
end

% =========================================================================
% Parameter and validation helpers
% =========================================================================
function params = make_params(cfg, sigma_X_sq, SNR_dB)
    params = struct();
    params.M          = cfg.M;
    params.P_avg      = cfg.P_avg;
    params.SNR_dB     = SNR_dB;
    params.sigma_X_sq = sigma_X_sq;
    params.sigma_n_sq = cfg.P_avg^2 / (10^(SNR_dB/10));
    params.R          = cfg.R;

    sig_t_sq = log(1 + sigma_X_sq);
    params.sig_t = sqrt(sig_t_sq);
    params.mu_t  = -0.5 * sig_t_sq;
end

function strictSummary = strict_validate_solution(x, p, params, x_pam, p_uni, strictConfigs)
    strictSummary = empty_strict_summary();

    ami_pam_val = AMI_functions.AMI_noCSI_validate(x_pam(:).', p_uni(:).', params);
    ami_standard = AMI_functions.AMI_noCSI_validate(x(:).', p(:).', params);

    strictResults = repmat(struct('name','','Ny',NaN,'yMin',NaN,'yMax',NaN, ...
        'ami',NaN,'gain',NaN,'runtime',NaN,'massErrMax',NaN,'masses',[]), numel(strictConfigs), 1);

    for k = 1:numel(strictConfigs)
        yg = build_strict_y_grid(params, x, strictConfigs(k));
        t0 = tic;
        [mi_strict, masses] = ami_validate_on_custom_grid(x, p, params, yg);
        rt = toc(t0);

        strictResults(k).name = strictConfigs(k).name;
        strictResults(k).Ny = numel(yg);
        strictResults(k).yMin = yg(1);
        strictResults(k).yMax = yg(end);
        strictResults(k).ami = mi_strict;
        strictResults(k).gain = mi_strict - ami_pam_val;
        strictResults(k).runtime = rt;
        strictResults(k).massErrMax = max(abs(masses - 1));
        strictResults(k).masses = masses;
    end

    strictAMIs = [strictResults.ami];
    allVals = [ami_standard, strictAMIs];

    strictSummary.ami_standard = ami_standard;
    strictSummary.strictResults = strictResults;
    strictSummary.ami_min = min(allVals);
    strictSummary.ami_max = max(allVals);
    strictSummary.ami_mean = mean(allVals);
    strictSummary.gain_min = strictSummary.ami_min - ami_pam_val;
    strictSummary.gain_max = strictSummary.ami_max - ami_pam_val;
    strictSummary.maxMassErr = max([strictResults.massErrMax]);
    strictSummary.errorMessage = '';
end

function configs = make_strict_configs()
    configs = struct([]);

    configs(1).name = 'strict-A';
    configs(1).tailProb = 1e-9;
    configs(1).coreSigma = 18;
    configs(1).dyCoreDiv = 28;
    configs(1).nMid = 1800;
    configs(1).nTail = 1800;

    configs(2).name = 'strict-B wider';
    configs(2).tailProb = 1e-11;
    configs(2).coreSigma = 22;
    configs(2).dyCoreDiv = 35;
    configs(2).nMid = 2600;
    configs(2).nTail = 2600;
end

function y_grid = build_strict_y_grid(params, x, cfg)
    sig = sqrt(params.sigma_n_sq);
    R = params.R;
    xM = max(abs(x));

    if params.sig_t > 1e-12
        h_hi = exp(params.mu_t + norminv(1 - cfg.tailProb) * params.sig_t);
    else
        h_hi = 1;
    end

    yAbs = R * xM * h_hi + 12*sig;
    yAbs = max(yAbs, 60*sig);

    yCore = cfg.coreSigma * sig;
    dyCore = max(sig / cfg.dyCoreDiv, 1e-5);
    core = -yCore : dyCore : yCore;

    startTail = max(yCore + dyCore, 2*dyCore);
    if startTail < yAbs
        posTail = exp(linspace(log(startTail), log(yAbs), cfg.nTail));
    else
        posTail = [];
    end

    yMidMax = min(0.25*yAbs, max(8*yCore, yCore + 20*sig));
    if startTail < yMidMax
        posMid = exp(linspace(log(startTail), log(yMidMax), cfg.nMid));
    else
        posMid = [];
    end

    y_pos = unique([posMid, posTail]);
    y_grid = sort(unique([-fliplr(y_pos), core, y_pos]));
end

function [mi_bits, masses] = ami_validate_on_custom_grid(x, p, params, y_grid)
    x = x(:).';
    p = p(:).';
    M = numel(x);

    Pyx = zeros(M, numel(y_grid));
    for i = 1:M
        Pyx(i,:) = calculate_Py_given_x(y_grid, x(i), params);
    end

    masses = zeros(M,1);
    for i = 1:M
        masses(i) = trapz(y_grid, Pyx(i,:));
    end

    mi_bits = AMI_functions.calculate_mi_quadrature(Pyx, p, y_grid);
end

% =========================================================================
% Random feasible initialization
% =========================================================================
function [x, p] = make_random_feasible_start(M, P_avg, minGap, p_min, xMaxUB, seed)
    rng(seed, 'twister');
    maxAttempts = 500;

    for attempt = 1:maxAttempts
        alphaChoices = [0.35, 0.7, 1.0, 2.0];
        alpha = alphaChoices(1 + mod(seed + attempt - 1, numel(alphaChoices)));
        v = rand(M,1).^(1/alpha);
        p = p_min + (1 - M*p_min) * v / sum(v);

        span = P_avg + (xMaxUB - P_avg) * rand();
        raw = sort(span * rand(M,1).^1.5, 'ascend');
        raw(1) = 0;
        raw(end) = max(raw(end), P_avg + minGap*(M-1));
        raw = sort(raw, 'ascend');

        cfgLocal = struct();
        cfgLocal.P_avg = P_avg;
        cfgLocal.SA = struct();
        cfgLocal.SA.minGap = minGap;
        cfgLocal.SA.enforce_sort = true;
        cfgLocal.SA.imdd_mode = true;
        cfgLocal.SA.enforce_power = true;
        cfgLocal.SA.powerConstraint = "mean";

        try
            x = project_constellation_weighted_power(raw, p, cfgLocal);
            x = x(:);
            p = project_probabilities_power(p(:), x, P_avg, p_min);

            ok = all(isfinite(x)) && all(isfinite(p)) && ...
                 all(x >= -1e-10) && all(p >= p_min - 1e-12) && ...
                 abs(sum(p)-1) < 1e-9 && abs(dot(p,x)-P_avg) < 1e-8 && ...
                 min(diff(x)) >= minGap - 1e-10 && max(x) <= xMaxUB + 1e-8;

            if ok
                x(abs(x) < 1e-12) = 0;
                return;
            end
        catch
        end
    end

    error('Could not generate a random feasible start after %d attempts.', maxAttempts);
end

% =========================================================================
% Summary tables and plots
% =========================================================================
function T = build_summary_table(results, M)
    n = numel(results);

    caseIndex = [results.caseIndex].';
    sigma_X_sq = [results.sigma_X_sq].';
    SNR_dB = [results.SNR_dB].';
    AMI_PAM_val = [results.ami_pam_val].';
    AMI_PAM_fast = [results.ami_pam_fast].';
    AMI_joint_fast = [results.ami_fast].';
    AMI_joint_val = [results.ami_val].';
    AMI_best_SA_val = nan(n,1);
    AMI_best_SA_fast = nan(n,1);
    SA_runtime_sec = nan(n,1);
    Gain_val = [results.gain_val].';
    bestStart = [results.bestStart].';
    runtime_sec = [results.runtime].';

    AMI_joint_strict_min = nan(n,1);
    AMI_joint_strict_mean = nan(n,1);
    Gain_strict_min = nan(n,1);
    strict_mass_err_max = nan(n,1);
    x_max = nan(n,1);
    min_p = nan(n,1);
    min_diff_x = nan(n,1);

    X = nan(n,M);
    P = nan(n,M);

    for i = 1:n
        if isfield(results(i).strict, 'ami_min')
            AMI_joint_strict_min(i) = results(i).strict.ami_min;
            AMI_joint_strict_mean(i) = results(i).strict.ami_mean;
            Gain_strict_min(i) = results(i).strict.gain_min;
            strict_mass_err_max(i) = results(i).strict.maxMassErr;
        end
        if ~isempty(results(i).allStarts)
            try
                sidx = results(i).bestStart;
                AMI_best_SA_val(i) = results(i).allStarts(sidx).ami_sa_val;
                AMI_best_SA_fast(i) = results(i).allStarts(sidx).ami_sa_fast;
                SA_runtime_sec(i) = results(i).allStarts(sidx).saRuntime;
            catch
            end
        end
        x = results(i).x_best(:).';
        p = results(i).p_best(:).';
        X(i,:) = x;
        P(i,:) = p;
        x_max(i) = max(x);
        min_p(i) = min(p);
        min_diff_x(i) = min(diff(x));
    end

    T = table(caseIndex, sigma_X_sq, SNR_dB, AMI_PAM_fast, AMI_PAM_val, ...
        AMI_best_SA_fast, AMI_best_SA_val, AMI_joint_fast, AMI_joint_val, ...
        AMI_joint_strict_min, AMI_joint_strict_mean, Gain_val, Gain_strict_min, ...
        bestStart, runtime_sec, SA_runtime_sec, x_max, min_p, min_diff_x, strict_mass_err_max);

    for k = 1:M
        T.(sprintf('x%d', k)) = X(:,k);
    end
    for k = 1:M
        T.(sprintf('p%d', k)) = P(:,k);
    end
end

function make_all_plots(T, turbulenceList, snrList, outDir)
    % Combined AMI plot.
    fig1 = figure('Color','w','Name','Joint sweep AMI combined');
    hold on; grid on; box on;
    for a = 1:numel(turbulenceList)
        sig = turbulenceList(a);
        idx = abs(T.sigma_X_sq - sig) < 1e-12;
        Ti = sortrows(T(idx,:), 'SNR_dB');
        plot(Ti.SNR_dB, Ti.AMI_joint_strict_min, 'o-', 'LineWidth', 1.4, ...
            'DisplayName', sprintf('Hybrid Joint-SA+fmincon GS/PS, sigma=%.1f', sig));
    end
    xlabel('SNR [dB]');
    ylabel('Strict validated AMI [bits/symbol]');
    title('Hybrid Joint-SA + Joint-fmincon AMI vs SNR');
    legend('Location','best');
    exportgraphics(fig1, fullfile(outDir, 'Joint_SA_fmincon_sweep_AMI_combined.png'), 'Resolution', 300);

    % Combined gain plot.
    fig2 = figure('Color','w','Name','Joint SA + fmincon sweep gain combined');
    hold on; grid on; box on;
    for a = 1:numel(turbulenceList)
        sig = turbulenceList(a);
        idx = abs(T.sigma_X_sq - sig) < 1e-12;
        Ti = sortrows(T(idx,:), 'SNR_dB');
        plot(Ti.SNR_dB, Ti.Gain_strict_min, 'o-', 'LineWidth', 1.4, ...
            'DisplayName', sprintf('sigma=%.1f', sig));
    end
    xlabel('SNR [dB]');
    ylabel('Conservative gain over Uniform PAM [bits/symbol]');
    title('Hybrid Joint-SA + Joint-fmincon conservative gain vs SNR');
    legend('Location','best');
    exportgraphics(fig2, fullfile(outDir, 'Joint_SA_fmincon_sweep_gain_combined.png'), 'Resolution', 300);

    % Per-turbulence plots with PAM baseline and optimized curve.
    for a = 1:numel(turbulenceList)
        sig = turbulenceList(a);
        idx = abs(T.sigma_X_sq - sig) < 1e-12;
        Ti = sortrows(T(idx,:), 'SNR_dB');
        safeSig = strrep(sprintf('%.1f', sig), '.', 'p');

        fA = figure('Color','w','Name',sprintf('AMI sigma %.1f', sig));
        plot(Ti.SNR_dB, Ti.AMI_PAM_val, 'o-', 'LineWidth', 1.3, 'DisplayName','Uniform PAM'); hold on;
        plot(Ti.SNR_dB, Ti.AMI_joint_strict_min, 'x-', 'LineWidth', 1.5, 'DisplayName','Hybrid Joint-SA+fmincon GS/PS');
        grid on; box on;
        xlabel('SNR [dB]');
        ylabel('AMI [bits/symbol]');
        title(sprintf('AMI vs SNR, sigma_X^2 = %.1f', sig));
        legend('Location','best');
        exportgraphics(fA, fullfile(outDir, sprintf('Joint_SA_fmincon_sweep_AMI_by_turbulence_sigma_%s.png', safeSig)), 'Resolution', 300);

        fG = figure('Color','w','Name',sprintf('Gain sigma %.1f', sig));
        plot(Ti.SNR_dB, Ti.Gain_strict_min, 'o-', 'LineWidth', 1.5);
        grid on; box on;
        xlabel('SNR [dB]');
        ylabel('Conservative gain [bits/symbol]');
        title(sprintf('Hybrid Joint-SA+fmincon GS/PS gain vs SNR, sigma_X^2 = %.1f', sig));
        exportgraphics(fG, fullfile(outDir, sprintf('Joint_SA_fmincon_sweep_gain_by_turbulence_sigma_%s.png', safeSig)), 'Resolution', 300);
    end
end

% =========================================================================
% Empty structs
% =========================================================================
function res = empty_case_result(M)
    res = struct();
    res.caseIndex = NaN;
    res.sigma_X_sq = NaN;
    res.SNR_dB = NaN;
    res.params = struct();
    res.x_pam = nan(M,1);
    res.p_uni = nan(M,1);
    res.ami_pam_fast = NaN;
    res.ami_pam_val = NaN;
    res.allStarts = repmat(empty_start_result(M), 0, 1);
    res.bestStart = NaN;
    res.x_best = nan(M,1);
    res.p_best = nan(M,1);
    res.ami_fast = NaN;
    res.ami_val = NaN;
    res.gain_val = NaN;
    res.strict = empty_strict_summary();
    res.runtime = NaN;
end

function sr = empty_start_result(M)
    sr = struct();
    sr.startIndex = NaN;
    sr.seed = NaN;
    sr.x0 = nan(M,1);
    sr.p0 = nan(M,1);
    sr.ami0_fast = NaN;
    sr.ami0_val = NaN;
    sr.x_sa = nan(M,1);
    sr.p_sa = nan(M,1);
    sr.ami_sa_fast = NaN;
    sr.ami_sa_val = NaN;
    sr.saRuntime = NaN;
    sr.saAcceptRate = NaN;
    sr.x_best = nan(M,1);
    sr.p_best = nan(M,1);
    sr.ami_fast = NaN;
    sr.ami_val = NaN;
    sr.exitflag = NaN;
    sr.ok = false;
    sr.errorMessage = '';
end

function st = empty_strict_summary()
    st = struct();
    st.ami_standard = NaN;
    st.strictResults = struct([]);
    st.ami_min = NaN;
    st.ami_max = NaN;
    st.ami_mean = NaN;
    st.gain_min = NaN;
    st.gain_max = NaN;
    st.maxMassErr = NaN;
    st.errorMessage = '';
end
