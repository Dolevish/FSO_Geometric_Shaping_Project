%% run_GS_M8_SNR15_hybrid_refinement.m
% Single-case geometric shaping simulation with Hybrid SA + Local Refinement:
% M = 8, SNR = 15 dB
%
% Goal:
%   1. Build uniform positive PAM baseline.
%   2. Compute baseline AMI.
%   3. Run multi-start Simulated Annealing for geometric shaping.
%   4. Refine the best SA candidates using local constrained optimization.
%   5. Report the constellation with the highest validated AMI.
%
% Output:
%   - Best GS constellation after refinement
%   - Validated AMI of Uniform PAM
%   - Validated AMI of Hybrid GS
%   - Validated AMI gain = AMI_GS - AMI_PAM

clear; clc; close all;
addpath(pwd);

fprintf('============================================================\n');
fprintf('   Hybrid Geometric Shaping: SA + Local Refinement\n');
fprintf('   Case: M=8, SNR=15 dB\n');
fprintf('============================================================\n\n');

%% 1. Channel and simulation parameters

params.R          = 1;      % Detector responsivity
params.M          = 8;      % Constellation size
params.P_avg      = 1;      % Average optical power
params.SNR_dB     = 15;     % SNR in dB

% Choose turbulence level:
% sigma_X_sq = sigma_R^2 in the paper.
% Try: 0, 0.1, 0.2, 0.3
params.sigma_X_sq = 0.1;

% According to the paper: SNR = P_avg^2 / sigma_n^2
params.sigma_n_sq = params.P_avg^2 / 10^(params.SNR_dB/10);

% Log-normal turbulence parameters:
% t = ln(h) ~ N(mu_t, sig_t^2)
sig_t_sq     = log(1 + params.sigma_X_sq);
params.sig_t = sqrt(sig_t_sq);
params.mu_t  = -0.5 * sig_t_sq;

fprintf('[1] Parameters\n');
fprintf('    M             = %d\n', params.M);
fprintf('    SNR           = %.1f dB\n', params.SNR_dB);
fprintf('    P_avg         = %.4f\n', params.P_avg);
fprintf('    sigma_X_sq    = %.4f\n', params.sigma_X_sq);
fprintf('    sigma_n_sq    = %.6e\n', params.sigma_n_sq);
fprintf('    mu_t          = %.6f\n', params.mu_t);
fprintf('    sig_t         = %.6f\n\n', params.sig_t);


%% 2. Baseline constellation: Uniform positive PAM

imdd_mode       = true;
powerConstraint = "mean";

x_pam = define_constellation( ...
    params.M, ...
    params.P_avg, ...
    imdd_mode, ...
    powerConstraint);

px = ones(params.M, 1) / params.M;

fprintf('[2] Uniform PAM baseline\n');
fprintf('    x_pam = ['); fprintf(' %.6f', x_pam); fprintf(' ]\n');
fprintf('    mean(x_pam) = %.12f\n\n', mean(x_pam));


%% 3. Build fixed y-grid for fast AMI evaluator

ghN_h = 40;

% Build a grid wide enough for all candidate constellations.
% Since projection keeps mean power fixed but max(x) can move,
% we use a conservative bound.
x_max_bound = max(x_pam) * 3.0;

y_grid = AMI_functions.build_noCSI_y_grid(params, x_max_bound);

fprintf('[3] Fast AMI y-grid\n');
fprintf('    y range = [%.4f, %.4f]\n', y_grid(1), y_grid(end));
fprintf('    Ny      = %d\n\n', numel(y_grid));


%% 4. Configure SA

cfg = struct();
cfg.params = params;
cfg.M      = params.M;
cfg.P_avg  = params.P_avg;

% Fast objective used inside simulated annealing
cfg.AMI_Evaluator = @(x) AMI_functions.AMI_noCSI_fast_grid( ...
    x(:).', px, params, ghN_h, y_grid);

cfg.SA = struct();

% Constraint/projection settings
cfg.SA.imdd_mode       = true;
cfg.SA.enforce_sort    = true;
cfg.SA.enforce_power   = true;
cfg.SA.powerConstraint = "mean";

% Numerical regularization
cfg.SA.minGap       = 0;
cfg.SA.projectIters = 3;

% SA search settings
cfg.SA.nStarts       = 10;
cfg.SA.maxIter       = 10000;
cfg.SA.itersPerTemp  = 150;

cfg.SA.T0            = 0.5;
cfg.SA.Tf            = 1e-3;
cfg.SA.nBlocks       = ceil(cfg.SA.maxIter / cfg.SA.itersPerTemp);
cfg.SA.coolingRate   = exp(log(cfg.SA.Tf / cfg.SA.T0) / cfg.SA.nBlocks);

cfg.SA.baseStd0      = 0.20;
cfg.SA.baseStdMin    = 1e-4;
cfg.SA.baseStdMax    = 0.80;
cfg.SA.baseStdGrow   = 1.10;
cfg.SA.baseStdShrink = 0.80;

cfg.SA.targetAccLo   = 0.20;
cfg.SA.targetAccHi   = 0.60;

cfg.SA.logEvery      = 500;

% Reproducibility.
% Use [] for a random seed each run.
cfg.SA.seedInit      = 12345;

% Parallel settings
cfg.SA.useParallel       = false;
cfg.SA.numWorkers        = [];
cfg.SA.closePoolWhenDone = false;

fprintf('[4] SA configuration\n');
fprintf('    nStarts      = %d\n', cfg.SA.nStarts);
fprintf('    maxIter      = %d\n', cfg.SA.maxIter);
fprintf('    minGap       = %.4f\n', cfg.SA.minGap);
fprintf('    seedInit     = %s\n\n', mat2str(cfg.SA.seedInit));


%% 5. Baseline AMI

fprintf('[5] Computing baseline AMI...\n');

ami_pam_fast = AMI_functions.AMI_noCSI_fast_grid( ...
    x_pam, px, params, ghN_h, y_grid);

ami_pam_val = AMI_functions.AMI_noCSI_validate( ...
    x_pam, px, params);

fprintf('    PAM AMI fast      = %.8f bits/symbol\n', ami_pam_fast);
fprintf('    PAM AMI validated = %.8f bits/symbol\n\n', ami_pam_val);


%% 6. Run multi-start simulated annealing

fprintf('[6] Running geometric shaping optimization with SA...\n\n');

[out, results] = sa_multistart(cfg, x_pam(:), "[M8-SNR15]");

x_sa_best = out.bestX(:).';
ami_sa_fast_best = out.bestMI;
ami_sa_val_best  = AMI_functions.AMI_noCSI_validate(x_sa_best, px, params);
gain_sa_val_best = ami_sa_val_best - ami_pam_val;

fprintf('\n[6] Best SA-only result\n');
fprintf('    best start       = %d\n', out.bestStart);
fprintf('    x_sa_best        = ['); fprintf(' %.6f', x_sa_best); fprintf(' ]\n');
fprintf('    SA AMI fast      = %.8f bits/symbol\n', ami_sa_fast_best);
fprintf('    SA AMI validated = %.8f bits/symbol\n', ami_sa_val_best);
fprintf('    SA gain validated= %.8f bits/symbol\n\n', gain_sa_val_best);


%% 7. Hybrid local refinement after SA
%
% The SA stage is good at finding a promising region of the search space,
% but it may not converge exactly to the best continuous constellation.
% This section takes the top SA candidates and refines each one using fmincon.
%
% The local problem is:
%       maximize AMI_fast(x)
% subject to:
%       x_i >= 0
%       mean(x) = P_avg
%       x_{i+1} - x_i >= minGap
%
% Since fmincon minimizes, the objective is -AMI_fast(x).

fprintf('[7] Running local constrained refinement on top SA candidates...\n');

useLocalRefinement = true;
localCandidates = struct([]);

if useLocalRefinement && exist('fmincon', 'file') == 2
    topK = min(5, numel(results));
    [~, idxSort] = sort([results.bestMI], 'descend');
    topIdx = idxSort(1:topK);

    fprintf('    Refining topK = %d SA candidates using fmincon...\n\n', topK);

    for kk = 1:topK
        kStart = topIdx(kk);
        x0_local = results(kStart).x_best(:);

        fprintf('    [Refine %d/%d] From SA start #%d | SA fast AMI = %.8f\n', ...
            kk, topK, kStart, results(kStart).bestMI);

        [x_ref, ami_ref_fast, exitflag, output] = local_refine_constellation_fmincon( ...
            x0_local, cfg, px, params, ghN_h, y_grid);

        x_ref = x_ref(:).';
        ami_ref_val = AMI_functions.AMI_noCSI_validate(x_ref, px, params);
        gain_ref_val = ami_ref_val - ami_pam_val;

        localCandidates(kk).sourceStart   = kStart;
        localCandidates(kk).x             = x_ref;
        localCandidates(kk).ami_fast      = ami_ref_fast;
        localCandidates(kk).ami_val       = ami_ref_val;
        localCandidates(kk).gain_val      = gain_ref_val;
        localCandidates(kk).exitflag      = exitflag;
        localCandidates(kk).output        = output;

        fprintf('        Refined x          = ['); fprintf(' %.6f', x_ref); fprintf(' ]\n');
        fprintf('        Refined AMI fast   = %.8f bits/symbol\n', ami_ref_fast);
        fprintf('        Refined AMI val    = %.8f bits/symbol\n', ami_ref_val);
        fprintf('        Refined gain val   = %.8f bits/symbol\n', gain_ref_val);
        fprintf('        exitflag           = %d\n\n', exitflag);
    end

    % Select final candidate by validated AMI, not only fast AMI.
    [~, bestLocalIdx] = max([localCandidates.ami_val]);

    x_gs       = localCandidates(bestLocalIdx).x;
    ami_gs_fast = localCandidates(bestLocalIdx).ami_fast;
    ami_gs_val  = localCandidates(bestLocalIdx).ami_val;

    fprintf('[7] Best refined candidate selected by validated AMI\n');
    fprintf('    source SA start  = %d\n', localCandidates(bestLocalIdx).sourceStart);
    fprintf('    best refined x   = ['); fprintf(' %.6f', x_gs); fprintf(' ]\n');
    fprintf('    best AMI fast    = %.8f bits/symbol\n', ami_gs_fast);
    fprintf('    best AMI val     = %.8f bits/symbol\n\n', ami_gs_val);

else
    fprintf('    fmincon was not found, or local refinement is disabled.\n');
    fprintf('    Falling back to the best SA-only constellation.\n\n');

    localCandidates = [];
    x_gs = AMI_functions.project_constellation_1D(x_sa_best(:), cfg);
    x_gs = x_gs(:).';
    ami_gs_fast = AMI_functions.AMI_noCSI_fast_grid(x_gs, px, params, ghN_h, y_grid);
    ami_gs_val  = AMI_functions.AMI_noCSI_validate(x_gs, px, params);
end


%% 8. Final validation and gain report

fprintf('[8] Final validation and gain report\n');

% Project once more just to be safe
x_gs = AMI_functions.project_constellation_1D(x_gs(:), cfg);
x_gs = x_gs(:).';

ami_gs_fast = AMI_functions.AMI_noCSI_fast_grid( ...
    x_gs, px, params, ghN_h, y_grid);

ami_gs_val = AMI_functions.AMI_noCSI_validate( ...
    x_gs, px, params);

gain_fast = ami_gs_fast - ami_pam_fast;
gain_val  = ami_gs_val  - ami_pam_val;

fprintf('    Final GS x        = ['); fprintf(' %.6f', x_gs); fprintf(' ]\n');
fprintf('    mean(x_gs)        = %.12f\n', mean(x_gs));
fprintf('    min(x_gs)         = %.12f\n', min(x_gs));
fprintf('    min diff(x_gs)    = %.12f\n\n', min(diff(x_gs)));

fprintf('    PAM AMI fast      = %.8f bits/symbol\n', ami_pam_fast);
fprintf('    PAM AMI validated = %.8f bits/symbol\n', ami_pam_val);
fprintf('    GS AMI fast       = %.8f bits/symbol\n', ami_gs_fast);
fprintf('    GS AMI validated  = %.8f bits/symbol\n', ami_gs_val);
fprintf('    AMI gain fast     = %.8f bits/symbol\n', gain_fast);
fprintf('    AMI gain validated= %.8f bits/symbol\n\n', gain_val);


%% 9. Save results

result = struct();
result.params          = params;
result.cfg             = cfg;
result.x_pam           = x_pam;
result.x_sa_best       = x_sa_best;
result.x_gs            = x_gs;
result.ami_pam_fast    = ami_pam_fast;
result.ami_pam_val     = ami_pam_val;
result.ami_sa_fast     = ami_sa_fast_best;
result.ami_sa_val      = ami_sa_val_best;
result.gain_sa_val     = gain_sa_val_best;
result.ami_gs_fast     = ami_gs_fast;
result.ami_gs_val      = ami_gs_val;
result.gain_fast       = gain_fast;
result.gain_val        = gain_val;
result.out             = out;
result.results         = results;
result.localCandidates = localCandidates;

save('result_GS_M8_SNR15_hybrid_refinement.mat', 'result');

fprintf('[9] Saved results to result_GS_M8_SNR15_hybrid_refinement.mat\n\n');


%% 10. Optional plot

figure;
stem(1:params.M, x_pam, 'o', 'LineWidth', 1.5); hold on;
stem(1:params.M, x_sa_best, 's', 'LineWidth', 1.5);
stem(1:params.M, x_gs, 'x', 'LineWidth', 1.5);
grid on; box on;
xlabel('Symbol index');
ylabel('Constellation level x');
title(sprintf('M=8, SNR=%g dB, \\sigma_R^2=%.2f', ...
    params.SNR_dB, params.sigma_X_sq));
legend('Uniform PAM', 'SA only', 'SA + local refinement', 'Location', 'best');

fprintf('Done.\n');


%% ========================================================================
% Local helper function: fmincon refinement
% ========================================================================

function [x_ref, ami_ref_fast, exitflag, output] = local_refine_constellation_fmincon( ...
    x0, cfg, px, params, ghN_h, y_grid)
%LOCAL_REFINE_CONSTELLATION_FMINCON
% Refines a feasible 1D IM/DD constellation using fmincon.
%
% Optimization problem:
%   maximize AMI_fast(x)
% subject to:
%   x_i >= 0
%   mean(x) = P_avg
%   x_{i+1} - x_i >= minGap
%
% Since fmincon minimizes, objective = -AMI_fast(x).

    M      = params.M;
    minGap = cfg.SA.minGap;

    % Make sure initial point is feasible according to the same projection
    % used by the SA algorithm.
    x0 = AMI_functions.project_constellation_1D(x0(:), cfg);
    x0 = x0(:);

    % Objective: negative AMI, because fmincon is a minimizer.
    objective = @(x) -AMI_functions.AMI_noCSI_fast_grid( ...
        x(:).', px, params, ghN_h, y_grid);

    % Equality constraint: mean(x) = P_avg.
    Aeq = ones(1, M) / M;
    beq = params.P_avg;

    % Inequality constraints for sorting + minimum gap:
    % x_{i+1} - x_i >= minGap
    % Equivalent form for fmincon: x_i - x_{i+1} <= -minGap
    A = zeros(M-1, M);
    b = -minGap * ones(M-1, 1);
    for i = 1:(M-1)
        A(i, i)   = 1;
        A(i, i+1) = -1;
    end

    % Non-negativity.
    lb = zeros(M, 1);
    ub = [];

    % fmincon options.
    options = optimoptions('fmincon', ...
        'Display', 'none', ...
        'Algorithm', 'sqp', ...
        'MaxIterations', 300, ...
        'MaxFunctionEvaluations', 5000, ...
        'OptimalityTolerance', 1e-8, ...
        'StepTolerance', 1e-10, ...
        'ConstraintTolerance', 1e-10);

    [x_ref, negAMI, exitflag, output] = fmincon( ...
        objective, x0, A, b, Aeq, beq, lb, ub, [], options);

    % Numerical cleanup: project again to keep exactly the same constraints as SA.
    x_ref = AMI_functions.project_constellation_1D(x_ref(:), cfg);
    x_ref = x_ref(:);

    % Recompute fast AMI after projection.
    ami_ref_fast = AMI_functions.AMI_noCSI_fast_grid( ...
        x_ref(:).', px, params, ghN_h, y_grid);

    % If fmincon returned a value before the final projection, keep the
    % recomputed value. This avoids inconsistencies after cleanup.
    %#ok<NASGU>
    negAMI = -ami_ref_fast;
end
