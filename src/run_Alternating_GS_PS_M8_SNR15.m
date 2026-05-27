%% run_Alternating_GS_PS_M8_SNR15.m
% Alternating Geometric Shaping + Probabilistic Shaping case study.
%
% Case:
%   M = 8, SNR = 15 dB, sigma_R^2 = 0.1, P_avg = 1
%
% Alternating scheme:
%   Start from Uniform PAM + uniform probabilities.
%   For round = 1..maxRounds:
%       1) GS step: optimize x while p is fixed.
%          Important: because p is no longer uniform, the power constraint is
%          weighted: sum_i p_i*x_i = P_avg.
%       2) PS step: optimize p while x is fixed.
%          Constraints: p_i>=p_min, sum_i p_i=1, sum_i p_i*x_i=P_avg.
%
% Output:
%   - x, p, AMI, gain after each GS and PS step
%   - final best alternating GS/PS solution
%   - plots and MAT file

clear; clc; close all;
addpath(pwd);

fprintf('============================================================\n');
fprintf('   Alternating GS/PS: M=8, SNR=15 dB\n');
fprintf('============================================================\n\n');

%% ========================================================================
%  1. Case parameters
% =========================================================================
M              = 8;
P_avg          = 1;
SNR_dB         = 15;
sigma_X_sq     = 0.1;
minGap_GS      = 0.005;
p_min          = 1e-8;
ghN_h          = 40;
maxRounds      = 5;
stopTolValAMI  = 1e-4;
runPamPSReference = true;

params = struct();
params.R          = 1;
params.M          = M;
params.P_avg      = P_avg;
params.SNR_dB     = SNR_dB;
params.sigma_X_sq = sigma_X_sq;
params.sigma_n_sq = P_avg^2 / 10^(SNR_dB/10);

sig_t_sq     = log(1 + sigma_X_sq);
params.sig_t = sqrt(sig_t_sq);
params.mu_t  = -0.5 * sig_t_sq;

fprintf('[1] Case parameters\n');
fprintf('    M              = %d\n', M);
fprintf('    SNR            = %.1f dB\n', SNR_dB);
fprintf('    P_avg          = %.6f\n', P_avg);
fprintf('    sigma_X_sq     = %.6f\n', sigma_X_sq);
fprintf('    sigma_n_sq     = %.6e\n', params.sigma_n_sq);
fprintf('    mu_t           = %.6f\n', params.mu_t);
fprintf('    sig_t          = %.6f\n', params.sig_t);
fprintf('    minGap_GS      = %.6f\n', minGap_GS);
fprintf('    p_min          = %.2e\n', p_min);
fprintf('    maxRounds      = %d\n', maxRounds);
fprintf('    stopTolValAMI  = %.2e\n\n', stopTolValAMI);

%% ========================================================================
%  2. Baseline and grid
% =========================================================================
x_pam = define_constellation(M, P_avg, true, "mean");
p_uniform = ones(M,1)/M;

% Conservative grid bound for the alternating run.
% Previous M=8/SNR=15 solutions usually stay below ~5, but we use 8 for safety.
x_max_bound = 8;
y_grid = AMI_functions.build_noCSI_y_grid(params, x_max_bound);

ami_pam_fast = AMI_functions.AMI_noCSI_fast_grid(x_pam, p_uniform, params, ghN_h, y_grid);
ami_pam_val  = AMI_functions.AMI_noCSI_validate(x_pam, p_uniform, params);

fprintf('[2] Uniform PAM baseline\n');
fprintf('    x_pam = ['); fprintf(' %.6f', x_pam); fprintf(' ]\n');
fprintf('    p_uni = ['); fprintf(' %.6f', p_uniform); fprintf(' ]\n');
fprintf('    mean(x_pam)        = %.12f\n', mean(x_pam));
fprintf('    dot(p_uni,x_pam)   = %.12f\n', dot(p_uniform, x_pam(:)));
fprintf('    AMI fast           = %.8f bits/symbol\n', ami_pam_fast);
fprintf('    AMI validated      = %.8f bits/symbol\n', ami_pam_val);
fprintf('    y_grid range       = [%.4f, %.4f], Ny=%d\n\n', ...
    y_grid(1), y_grid(end), numel(y_grid));

%% ========================================================================
%  3. Shared config
% =========================================================================
cfg = struct();
cfg.params = params;
cfg.M      = M;
cfg.P_avg  = P_avg;

% GS evaluator: x varies, p is fixed and supplied by the alternating loop.
cfg.GS_Evaluator = @(x, p) AMI_functions.AMI_noCSI_fast_grid( ...
    x(:).', p(:).', params, ghN_h, y_grid);

% PS evaluator: x is fixed, p varies.
cfg.PS_Evaluator = @(x_fixed, p) AMI_functions.AMI_noCSI_fast_grid( ...
    x_fixed(:).', p(:).', params, ghN_h, y_grid);

%% ========================================================================
%  4. Configure weighted GS SA
% =========================================================================
cfg.SA = struct();
cfg.SA.imdd_mode       = true;
cfg.SA.enforce_sort    = true;
cfg.SA.enforce_power   = true;
cfg.SA.powerConstraint = "weighted_mean"; % documentation only; custom projection is used
cfg.SA.minGap          = minGap_GS;
cfg.SA.projectIters    = 3;

cfg.SA.nStarts       = 12;
cfg.SA.maxIter       = 8000;
cfg.SA.itersPerTemp  = 100;
cfg.SA.T0            = 0.4;
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
cfg.SA.logEvery      = 0;
cfg.SA.seedInit      = 12345;
cfg.SA.useParallel   = true;
cfg.SA.numWorkers    = [];
cfg.SA.closePoolWhenDone = false;

%% ========================================================================
%  5. Configure PS SA
% =========================================================================
cfg.PS = struct();
cfg.PS.p_min         = p_min;
cfg.PS.nStarts       = 12;
cfg.PS.maxIter       = 8000;
cfg.PS.itersPerTemp  = 100;

cfg.PS.T0            = 0.05;
cfg.PS.Tf            = 1e-5;
cfg.PS.nBlocks       = ceil(cfg.PS.maxIter / cfg.PS.itersPerTemp);
cfg.PS.coolingRate   = exp(log(cfg.PS.Tf / cfg.PS.T0) / cfg.PS.nBlocks);

cfg.PS.baseStd0      = 0.05;
cfg.PS.baseStdMin    = 1e-5;
cfg.PS.baseStdMax    = 0.50;
cfg.PS.baseStdGrow   = 1.10;
cfg.PS.baseStdShrink = 0.80;
cfg.PS.targetAccLo   = 0.20;
cfg.PS.targetAccHi   = 0.60;
cfg.PS.initNoiseStd  = 0.05;
cfg.PS.logEvery      = 0;
cfg.PS.seedInit      = 54321;
cfg.PS.useParallel   = true;
cfg.PS.numWorkers    = [];
cfg.PS.closePoolWhenDone = false;

fprintf('[3] Optimization settings\n');
fprintf('    Weighted GS: nStarts=%d, maxIter=%d, minGap=%.6f\n', ...
    cfg.SA.nStarts, cfg.SA.maxIter, cfg.SA.minGap);
fprintf('    PS:          nStarts=%d, maxIter=%d, p_min=%.2e\n\n', ...
    cfg.PS.nStarts, cfg.PS.maxIter, cfg.PS.p_min);

%% ========================================================================
%  6. Optional reference: Uniform PAM + PS
% =========================================================================
ref = struct();
if runPamPSReference
    fprintf('============================================================\n');
    fprintf(' [Reference] Uniform PAM + PS\n');
    fprintf('============================================================\n');

    [out_ps_pam, results_ps_pam] = ps_multistart(cfg, x_pam(:), p_uniform(:), '[PAM-PS]');
    p_pam_ps = project_probabilities_power(out_ps_pam.bestP(:), x_pam(:), P_avg, p_min);

    ref.p_pam_ps = p_pam_ps;
    ref.ami_pam_ps_fast = cfg.PS_Evaluator(x_pam(:), p_pam_ps);
    ref.ami_pam_ps_val  = AMI_functions.AMI_noCSI_validate(x_pam, p_pam_ps, params);
    ref.results_ps_pam  = results_ps_pam;

    print_result_block('Uniform PAM + PS', x_pam(:), p_pam_ps, ...
        ref.ami_pam_ps_fast, ref.ami_pam_ps_val, ref.ami_pam_ps_val - ami_pam_val, P_avg);
end

%% ========================================================================
%  7. Alternating GS/PS loop
% =========================================================================
fprintf('============================================================\n');
fprintf(' [Alternating GS/PS]\n');
fprintf('============================================================\n');

x_current = x_pam(:);
p_current = p_uniform(:);
current_val = ami_pam_val;

history = struct([]);
roundCounter = 0;

% Store round 0.
roundCounter = roundCounter + 1;
history(roundCounter).round = 0;
history(roundCounter).stage = "Initial Uniform PAM";
history(roundCounter).x = x_current(:);
history(roundCounter).p = p_current(:);
history(roundCounter).amiFast = ami_pam_fast;
history(roundCounter).amiVal = ami_pam_val;
history(roundCounter).gainVal = 0;
history(roundCounter).source = "baseline";

best = history(roundCounter);

for r = 1:maxRounds
    fprintf('\n------------------------------------------------------------\n');
    fprintf(' Alternating round %d/%d\n', r, maxRounds);
    fprintf('------------------------------------------------------------\n');

    %% ----- GS step with p fixed -----
    fprintf('\n[A] Weighted GS step: optimize x with p fixed\n');
    fprintf('    current p = ['); fprintf(' %.6e', p_current); fprintf(' ]\n');

    cfg.SA.seedInit = uint32(12345 + 1000*r);
    labelGS = sprintf('[ALT%02d-GS]', r);
    [out_wgs, results_wgs] = gs_weighted_multistart(cfg, x_current(:), p_current(:), labelGS);

    x_gs_sa = project_constellation_weighted_power(out_wgs.bestX(:), p_current(:), cfg);
    ami_gs_sa_fast = cfg.GS_Evaluator(x_gs_sa, p_current);
    ami_gs_sa_val  = AMI_functions.AMI_noCSI_validate(x_gs_sa, p_current, params);

    fprintf('\n    Best weighted-GS SA-only result\n');
    fprintf('    x_gs_sa = ['); fprintf(' %.6f', x_gs_sa); fprintf(' ]\n');
    fprintf('    dot(p,x) = %.12f\n', dot(p_current(:), x_gs_sa(:)));
    fprintf('    AMI fast = %.8f | AMI val = %.8f | gain = %.8f\n', ...
        ami_gs_sa_fast, ami_gs_sa_val, ami_gs_sa_val - ami_pam_val);

    refineTopK = 5;
    [x_gs, gs_ref_info] = refine_weighted_gs_top_candidates( ...
        cfg, results_wgs, p_current, params, ghN_h, y_grid, ami_pam_val, refineTopK);

    ami_gs_fast = cfg.GS_Evaluator(x_gs, p_current);
    ami_gs_val  = AMI_functions.AMI_noCSI_validate(x_gs, p_current, params);

    % Monotonic-by-validation guard: do not accept a GS step that lowers validated AMI.
    if ami_gs_val + 1e-10 < current_val
        fprintf('    Weighted-GS validation did not improve. Keeping previous x.\n');
        x_gs = x_current;
        ami_gs_fast = cfg.GS_Evaluator(x_gs, p_current);
        ami_gs_val  = current_val;
        gsAccepted = false;
    else
        gsAccepted = true;
    end

    print_result_block(sprintf('Round %d after weighted GS', r), x_gs(:), p_current(:), ...
        ami_gs_fast, ami_gs_val, ami_gs_val - ami_pam_val, P_avg);

    roundCounter = roundCounter + 1;
    history(roundCounter).round = r;
    history(roundCounter).stage = "After weighted GS";
    history(roundCounter).x = x_gs(:);
    history(roundCounter).p = p_current(:);
    history(roundCounter).amiFast = ami_gs_fast;
    history(roundCounter).amiVal = ami_gs_val;
    history(roundCounter).gainVal = ami_gs_val - ami_pam_val;
    history(roundCounter).source = string(sprintf('GS accepted=%d', gsAccepted));
    history(roundCounter).results_wgs = results_wgs;
    history(roundCounter).gs_ref_info = gs_ref_info;

    if history(roundCounter).amiVal > best.amiVal
        best = history(roundCounter);
    end

    %% ----- PS step with x fixed -----
    fprintf('\n[B] PS step: optimize p with x fixed\n');
    cfg.PS.seedInit = uint32(54321 + 1000*r);
    labelPS = sprintf('[ALT%02d-PS]', r);
    [out_ps, results_ps] = ps_multistart(cfg, x_gs(:), p_current(:), labelPS);

    p_ps = project_probabilities_power(out_ps.bestP(:), x_gs(:), P_avg, p_min);
    ami_ps_fast = cfg.PS_Evaluator(x_gs(:), p_ps);
    ami_ps_val  = AMI_functions.AMI_noCSI_validate(x_gs, p_ps, params);

    % Monotonic-by-validation guard.
    if ami_ps_val + 1e-10 < ami_gs_val
        fprintf('    PS validation did not improve. Keeping previous p.\n');
        p_ps = p_current;
        ami_ps_fast = cfg.PS_Evaluator(x_gs(:), p_ps);
        ami_ps_val = ami_gs_val;
        psAccepted = false;
    else
        psAccepted = true;
    end

    print_result_block(sprintf('Round %d after PS', r), x_gs(:), p_ps(:), ...
        ami_ps_fast, ami_ps_val, ami_ps_val - ami_pam_val, P_avg);

    roundCounter = roundCounter + 1;
    history(roundCounter).round = r;
    history(roundCounter).stage = "After PS";
    history(roundCounter).x = x_gs(:);
    history(roundCounter).p = p_ps(:);
    history(roundCounter).amiFast = ami_ps_fast;
    history(roundCounter).amiVal = ami_ps_val;
    history(roundCounter).gainVal = ami_ps_val - ami_pam_val;
    history(roundCounter).source = string(sprintf('PS accepted=%d', psAccepted));
    history(roundCounter).results_ps = results_ps;

    if history(roundCounter).amiVal > best.amiVal
        best = history(roundCounter);
    end

    improvement = ami_ps_val - current_val;
    fprintf('\n[C] Round %d improvement over previous alternating state: %.8f bits/symbol\n', ...
        r, improvement);

    x_current = x_gs(:);
    p_current = p_ps(:);
    current_val = ami_ps_val;

    if improvement >= 0 && improvement < stopTolValAMI
        fprintf('    Stopping: improvement %.3e < stopTol %.3e\n', improvement, stopTolValAMI);
        break;
    end
end

%% ========================================================================
%  8. Final summary
% =========================================================================
fprintf('\n============================================================\n');
fprintf(' Final alternating GS/PS summary\n');
fprintf('============================================================\n');
fprintf('---------------------------------------------------------------------------------------------\n');
fprintf(' idx | round | stage                  | AMI fast   | AMI validated | Gain validated\n');
fprintf('---------------------------------------------------------------------------------------------\n');
for i = 1:numel(history)
    fprintf(' %3d | %5d | %-22s | %10.8f | %13.8f | %14.8f\n', ...
        i, history(i).round, history(i).stage, history(i).amiFast, history(i).amiVal, history(i).gainVal);
end
fprintf('---------------------------------------------------------------------------------------------\n\n');

fprintf('Best alternating solution selected by validated AMI:\n');
fprintf('    stage      = %s\n', best.stage);
fprintf('    round      = %d\n', best.round);
fprintf('    AMI val    = %.8f bits/symbol\n', best.amiVal);
fprintf('    gain val   = %.8f bits/symbol\n', best.gainVal);
fprintf('    x_best     = ['); fprintf(' %.6f', best.x(:)); fprintf(' ]\n');
fprintf('    p_best     = ['); fprintf(' %.6e', best.p(:)); fprintf(' ]\n');
fprintf('    sum(p)     = %.12f\n', sum(best.p));
fprintf('    dot(p,x)   = %.12f\n', dot(best.p(:), best.x(:)));
fprintf('    min diff x = %.12e\n', min(diff(best.x(:))));
fprintf('    min(p)     = %.3e\n\n', min(best.p));

%% ========================================================================
%  9. Save and plot
% =========================================================================
result = struct();
result.params = params;
result.cfg = cfg;
result.M = M;
result.P_avg = P_avg;
result.SNR_dB = SNR_dB;
result.sigma_X_sq = sigma_X_sq;
result.minGap_GS = minGap_GS;
result.p_min = p_min;
result.ghN_h = ghN_h;
result.y_grid = y_grid;

result.x_pam = x_pam(:);
result.p_uniform = p_uniform(:);
result.ami_pam_fast = ami_pam_fast;
result.ami_pam_val = ami_pam_val;
result.ref = ref;
result.history = history;
result.best = best;

save('result_Alternating_GS_PS_M8_SNR15.mat', 'result');
fprintf('[Save] Saved result_Alternating_GS_PS_M8_SNR15.mat\n');

stages = string({history.stage});
amiVals = [history.amiVal];
gainVals = [history.gainVal];

fig1 = figure('Name','Alternating GS/PS AMI Progress','Color','w');
plot(0:numel(amiVals)-1, amiVals, '-o', 'LineWidth', 1.7, 'MarkerSize', 6);
grid on; box on;
xlabel('Alternating step index');
ylabel('Validated AMI (bits/symbol)');
title('Alternating GS/PS progress: M=8, SNR=15 dB');
xticks(0:numel(amiVals)-1);
xticklabels(stages);
xtickangle(30);
exportgraphics(fig1, 'Alternating_GS_PS_AMI_progress_M8_SNR15.png', 'Resolution', 200);

fig2 = figure('Name','Alternating GS/PS Gain Progress','Color','w');
plot(0:numel(gainVals)-1, gainVals, '-o', 'LineWidth', 1.7, 'MarkerSize', 6);
grid on; box on;
xlabel('Alternating step index');
ylabel('Validated gain over Uniform PAM (bits/symbol)');
title('Alternating GS/PS gain progress: M=8, SNR=15 dB');
xticks(0:numel(gainVals)-1);
xticklabels(stages);
xtickangle(30);
exportgraphics(fig2, 'Alternating_GS_PS_gain_progress_M8_SNR15.png', 'Resolution', 200);

fig3 = figure('Name','Best Alternating x and p','Color','w');
yyaxis left;
stem(1:M, best.x(:), 'o', 'LineWidth', 1.5); grid on; box on;
ylabel('Constellation level x_i');
yyaxis right;
stem(1:M, best.p(:), 'x', 'LineWidth', 1.5);
ylabel('Probability p_i');
xlabel('Symbol index');
title('Best Alternating GS/PS solution');
legend('x_i', 'p_i', 'Location', 'best');
exportgraphics(fig3, 'Alternating_GS_PS_best_x_p_M8_SNR15.png', 'Resolution', 200);

fprintf('[Plot] Saved Alternating_GS_PS_AMI_progress_M8_SNR15.png\n');
fprintf('[Plot] Saved Alternating_GS_PS_gain_progress_M8_SNR15.png\n');
fprintf('[Plot] Saved Alternating_GS_PS_best_x_p_M8_SNR15.png\n');
fprintf('Done.\n');


%% ========================================================================
%  Local helper functions
% =========================================================================
function print_result_block(name, x, p, ami_fast, ami_val, gain_val, P_avg)
    fprintf('\n    %s\n', name);
    fprintf('    x = ['); fprintf(' %.6f', x(:)); fprintf(' ]\n');
    fprintf('    p = ['); fprintf(' %.6e', p(:)); fprintf(' ]\n');
    fprintf('    sum(p)      = %.12f\n', sum(p));
    fprintf('    dot(p,x)    = %.12f\n', dot(p(:), x(:)));
    fprintf('    target power= %.12f\n', P_avg);
    fprintf('    min(p)      = %.3e\n', min(p));
    if numel(x) > 1
        fprintf('    min diff(x) = %.3e\n', min(diff(x(:))));
    end
    fprintf('    AMI fast    = %.8f bits/symbol\n', ami_fast);
    fprintf('    AMI val     = %.8f bits/symbol\n', ami_val);
    fprintf('    gain val    = %.8f bits/symbol\n\n', gain_val);
end


function [x_best_final, info] = refine_weighted_gs_top_candidates(cfg, results_gs, p_fixed, params, ghN_h, y_grid, ami_pam_val, topK)
    info = struct();
    info.usedFmincon = false;
    info.candidates = [];

    allMI = [results_gs.bestMI];
    [~, idxSort] = sort(allMI, 'descend');
    topK = min(topK, numel(idxSort));

    % Start with best SA candidate by validated AMI.
    bestVal = -Inf;
    x_best_final = results_gs(idxSort(1)).x_best(:);
    p_fixed = p_fixed(:);

    for t = 1:topK
        idx = idxSort(t);
        x0 = project_constellation_weighted_power(results_gs(idx).x_best(:), p_fixed, cfg);
        val0 = AMI_functions.AMI_noCSI_validate(x0, p_fixed, params);
        if val0 > bestVal
            bestVal = val0;
            x_best_final = x0(:);
        end
    end

    if exist('fmincon', 'file') ~= 2
        fprintf('    fmincon unavailable. Keeping best weighted-GS SA candidate by validated AMI.\n');
        return;
    end

    fprintf('    Refining topK=%d weighted-GS candidates using fmincon...\n', topK);
    info.usedFmincon = true;

    M = cfg.M;
    minGap = cfg.SA.minGap;

    Aeq = p_fixed(:).';
    beq = cfg.P_avg;
    lb = zeros(M,1);
    ub = [];

    A = zeros(M-1, M);
    b = -minGap * ones(M-1,1);
    for i = 1:M-1
        A(i,i)   = 1;
        A(i,i+1) = -1;
    end

    objective = @(x) -AMI_functions.AMI_noCSI_fast_grid(x(:).', p_fixed(:).', params, ghN_h, y_grid);

    options = optimoptions('fmincon', ...
        'Display','none', ...
        'Algorithm','sqp', ...
        'MaxIterations',300, ...
        'MaxFunctionEvaluations',5000, ...
        'OptimalityTolerance',1e-9, ...
        'ConstraintTolerance',1e-10, ...
        'StepTolerance',1e-10);

    cand = repmat(struct('sourceStart',NaN,'x',[],'amiFast',NaN,'amiVal',NaN, ...
        'gainVal',NaN,'exitflag',NaN), topK, 1);

    for t = 1:topK
        idx = idxSort(t);
        x0 = project_constellation_weighted_power(results_gs(idx).x_best(:), p_fixed, cfg);

        try
            [x_ref, negMI, exitflag] = fmincon(objective, x0(:), A, b, Aeq, beq, lb, ub, [], options);
            x_ref = project_constellation_weighted_power(x_ref(:), p_fixed, cfg);
            amiFast = -negMI;
            amiVal = AMI_functions.AMI_noCSI_validate(x_ref, p_fixed, params);
        catch ME
            fprintf('        Weighted-GS refinement failed for candidate %d: %s\n', t, ME.message);
            x_ref = x0(:);
            exitflag = NaN;
            amiFast = AMI_functions.AMI_noCSI_fast_grid(x_ref(:).', p_fixed(:).', params, ghN_h, y_grid);
            amiVal = AMI_functions.AMI_noCSI_validate(x_ref, p_fixed, params);
        end

        cand(t).sourceStart = results_gs(idx).startIndex;
        cand(t).x = x_ref(:);
        cand(t).amiFast = amiFast;
        cand(t).amiVal = amiVal;
        cand(t).gainVal = amiVal - ami_pam_val;
        cand(t).exitflag = exitflag;

        fprintf('        [Weighted GS refine %d/%d] start #%d | AMI val=%.8f | gain=%.8f | exitflag=%g\n', ...
            t, topK, cand(t).sourceStart, cand(t).amiVal, cand(t).gainVal, cand(t).exitflag);

        if amiVal > bestVal
            bestVal = amiVal;
            x_best_final = x_ref(:);
        end
    end

    info.candidates = cand;
end
