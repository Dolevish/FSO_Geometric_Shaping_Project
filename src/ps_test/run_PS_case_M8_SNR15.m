%% run_PS_case_M8_SNR15.m
% Case study: compare GS and PS for M=8, SNR=15 dB.
%
% Experiments:
%   1. Uniform PAM
%   2. Uniform PAM + PS
%   3. GS with positive minGap
%   4. GS with positive minGap + PS
%
% Notes:
%   - GS optimizes constellation levels x with uniform probabilities.
%   - PS optimizes probabilities p for a fixed constellation x.
%   - PS enforces:
%         p_i >= p_min
%         sum_i p_i = 1
%         sum_i p_i x_i = P_avg

clear; clc; close all;
addpath(pwd);

fprintf('============================================================\n');
fprintf('   PS Case Study: M=8, SNR=15 dB, Uniform PAM / GS / PS\n');
fprintf('============================================================\n\n');

%% ========================================================================
%  1. Fixed case parameters
% =========================================================================
M              = 8;
P_avg          = 1;
SNR_dB         = 15;
sigma_X_sq     = 0.1;
minGap_GS      = 0.005;
p_min          = 1e-8;
ghN_h          = 40;

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
fprintf('    p_min          = %.2e\n\n', p_min);

%% ========================================================================
%  2. Baseline constellation and grid
% =========================================================================
x_pam = define_constellation(M, P_avg, true, "mean");
p_uniform = ones(M,1)/M;

x_max_bound = max(x_pam) * 3.0;
y_grid = AMI_functions.build_noCSI_y_grid(params, x_max_bound);

fprintf('[2] Uniform PAM baseline\n');
fprintf('    x_pam = ['); fprintf(' %.6f', x_pam); fprintf(' ]\n');
fprintf('    p_uni = ['); fprintf(' %.6f', p_uniform); fprintf(' ]\n');
fprintf('    mean(x_pam)        = %.12f\n', mean(x_pam));
fprintf('    dot(p_uni,x_pam)   = %.12f\n', dot(p_uniform, x_pam(:)));
fprintf('    y_grid range       = [%.4f, %.4f], Ny=%d\n\n', ...
    y_grid(1), y_grid(end), numel(y_grid));

%% ========================================================================
%  3. Build shared config
% =========================================================================
cfg = struct();
cfg.params = params;
cfg.M      = M;
cfg.P_avg  = P_avg;

% Evaluator for GS: x varies, probabilities are uniform.
cfg.AMI_Evaluator = @(x) AMI_functions.AMI_noCSI_fast_grid( ...
    x(:).', p_uniform(:).', params, ghN_h, y_grid);

% Evaluator for PS: x is fixed, p varies.
cfg.PS_Evaluator = @(x_fixed, p) AMI_functions.AMI_noCSI_fast_grid( ...
    x_fixed(:).', p(:).', params, ghN_h, y_grid);

%% ========================================================================
%  4. Configure GS SA
% =========================================================================
cfg.SA = struct();
cfg.SA.imdd_mode       = true;
cfg.SA.enforce_sort    = true;
cfg.SA.enforce_power   = true;
cfg.SA.powerConstraint = "mean";
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
fprintf('    GS: nStarts=%d, maxIter=%d, minGap=%.6f\n', ...
    cfg.SA.nStarts, cfg.SA.maxIter, cfg.SA.minGap);
fprintf('    PS: nStarts=%d, maxIter=%d, p_min=%.2e\n\n', ...
    cfg.PS.nStarts, cfg.PS.maxIter, cfg.PS.p_min);

%% ========================================================================
%  6. Experiment 1: Uniform PAM
% =========================================================================
fprintf('============================================================\n');
fprintf(' [Experiment 1] Uniform PAM\n');
fprintf('============================================================\n');

ami_pam_fast = AMI_functions.AMI_noCSI_fast_grid(x_pam, p_uniform, params, ghN_h, y_grid);
ami_pam_val  = AMI_functions.AMI_noCSI_validate(x_pam, p_uniform, params);

fprintf('    AMI fast      = %.8f bits/symbol\n', ami_pam_fast);
fprintf('    AMI validated = %.8f bits/symbol\n\n', ami_pam_val);

%% ========================================================================
%  7. Experiment 2: Uniform PAM + PS
% =========================================================================
fprintf('============================================================\n');
fprintf(' [Experiment 2] Uniform PAM + PS\n');
fprintf('============================================================\n');

[out_ps_pam, results_ps_pam] = ps_multistart(cfg, x_pam(:), p_uniform(:), '[PAM-PS]');

p_pam_ps = out_ps_pam.bestP(:);
p_pam_ps = project_probabilities_power(p_pam_ps, x_pam(:), P_avg, p_min);

ami_pam_ps_fast = cfg.PS_Evaluator(x_pam(:), p_pam_ps);
ami_pam_ps_val  = AMI_functions.AMI_noCSI_validate(x_pam, p_pam_ps, params);
gain_pam_ps_val = ami_pam_ps_val - ami_pam_val;

print_result_block('Uniform PAM + PS', x_pam(:), p_pam_ps, ami_pam_ps_fast, ami_pam_ps_val, gain_pam_ps_val, P_avg);

%% ========================================================================
%  8. Experiment 3: GS with positive minGap
% =========================================================================
fprintf('============================================================\n');
fprintf(' [Experiment 3] GS with positive minGap\n');
fprintf('============================================================\n');

[out_gs, results_gs] = sa_multistart(cfg, x_pam(:), '[GS]');

x_gs_sa = out_gs.bestX(:);
ami_gs_sa_fast = cfg.AMI_Evaluator(x_gs_sa);
ami_gs_sa_val  = AMI_functions.AMI_noCSI_validate(x_gs_sa, p_uniform, params);

fprintf('\n    Best SA-only GS:\n');
fprintf('    x_gs_sa = ['); fprintf(' %.6f', x_gs_sa); fprintf(' ]\n');
fprintf('    AMI fast      = %.8f\n', ami_gs_sa_fast);
fprintf('    AMI validated = %.8f\n', ami_gs_sa_val);
fprintf('    gain val      = %.8f\n\n', ami_gs_sa_val - ami_pam_val);

% Optional but recommended: local constrained refinement for GS.
refineTopK = 5;
[x_gs, gs_ref_info] = refine_gs_top_candidates(cfg, results_gs, p_uniform, params, ghN_h, y_grid, ami_pam_val, refineTopK);

ami_gs_fast = cfg.AMI_Evaluator(x_gs);
ami_gs_val  = AMI_functions.AMI_noCSI_validate(x_gs, p_uniform, params);
gain_gs_val = ami_gs_val - ami_pam_val;

print_result_block('GS minGap>0', x_gs(:), p_uniform(:), ami_gs_fast, ami_gs_val, gain_gs_val, P_avg);

%% ========================================================================
%  9. Experiment 4: GS with positive minGap + PS
% =========================================================================
fprintf('============================================================\n');
fprintf(' [Experiment 4] GS with positive minGap + PS\n');
fprintf('============================================================\n');

[out_ps_gs, results_ps_gs] = ps_multistart(cfg, x_gs(:), p_uniform(:), '[GS-PS]');

p_gs_ps = out_ps_gs.bestP(:);
p_gs_ps = project_probabilities_power(p_gs_ps, x_gs(:), P_avg, p_min);

ami_gs_ps_fast = cfg.PS_Evaluator(x_gs(:), p_gs_ps);
ami_gs_ps_val  = AMI_functions.AMI_noCSI_validate(x_gs, p_gs_ps, params);
gain_gs_ps_val = ami_gs_ps_val - ami_pam_val;

print_result_block('GS minGap>0 + PS', x_gs(:), p_gs_ps, ami_gs_ps_fast, ami_gs_ps_val, gain_gs_ps_val, P_avg);

%% ========================================================================
%  10. Final comparison table
% =========================================================================
methods = [
    "Uniform PAM"
    "Uniform PAM + PS"
    "GS minGap>0"
    "GS minGap>0 + PS"
];

AMI_fast = [
    ami_pam_fast
    ami_pam_ps_fast
    ami_gs_fast
    ami_gs_ps_fast
];

AMI_val = [
    ami_pam_val
    ami_pam_ps_val
    ami_gs_val
    ami_gs_ps_val
];

Gain_val = AMI_val - ami_pam_val;

fprintf('\n============================================================\n');
fprintf(' Final comparison: M=8, SNR=15 dB, sigma_X_sq=0.1\n');
fprintf('============================================================\n');
fprintf('--------------------------------------------------------------------------------\n');
fprintf(' idx | Method              | AMI fast   | AMI validated | Gain validated\n');
fprintf('--------------------------------------------------------------------------------\n');
for i = 1:numel(methods)
    fprintf(' %3d | %-19s | %10.8f | %13.8f | %14.8f\n', ...
        i, methods(i), AMI_fast(i), AMI_val(i), Gain_val(i));
end
fprintf('--------------------------------------------------------------------------------\n\n');

fprintf('Final x/p values:\n');
fprintf('1) Uniform PAM\n');
print_x_p(x_pam(:), p_uniform(:), P_avg);
fprintf('2) Uniform PAM + PS\n');
print_x_p(x_pam(:), p_pam_ps(:), P_avg);
fprintf('3) GS minGap>0\n');
print_x_p(x_gs(:), p_uniform(:), P_avg);
fprintf('4) GS minGap>0 + PS\n');
print_x_p(x_gs(:), p_gs_ps(:), P_avg);

%% ========================================================================
%  11. Save and plot
% =========================================================================
result = struct();
result.params = params;
result.cfg = cfg;
result.minGap_GS = minGap_GS;
result.p_min = p_min;
result.ghN_h = ghN_h;
result.y_grid = y_grid;

result.x_pam = x_pam(:);
result.p_uniform = p_uniform(:);
result.ami_pam_fast = ami_pam_fast;
result.ami_pam_val = ami_pam_val;

result.p_pam_ps = p_pam_ps(:);
result.ami_pam_ps_fast = ami_pam_ps_fast;
result.ami_pam_ps_val = ami_pam_ps_val;
result.results_ps_pam = results_ps_pam;

result.x_gs = x_gs(:);
result.ami_gs_fast = ami_gs_fast;
result.ami_gs_val = ami_gs_val;
result.results_gs = results_gs;
result.gs_ref_info = gs_ref_info;

result.p_gs_ps = p_gs_ps(:);
result.ami_gs_ps_fast = ami_gs_ps_fast;
result.ami_gs_ps_val = ami_gs_ps_val;
result.results_ps_gs = results_ps_gs;

result.methods = methods;
result.AMI_fast = AMI_fast;
result.AMI_val = AMI_val;
result.Gain_val = Gain_val;

save('result_PS_case_M8_SNR15.mat', 'result');
fprintf('[Save] Saved result_PS_case_M8_SNR15.mat\n');

fig = figure('Name','PS Case Study AMI Comparison','Color','w');
bar(AMI_val);
grid on; box on;
set(gca, 'XTick', 1:numel(methods), 'XTickLabel', methods, 'XTickLabelRotation', 25);
ylabel('Validated AMI (bits/symbol)');
title('M=8, SNR=15 dB: Uniform PAM / PS / GS / GS+PS');
exportgraphics(fig, 'PS_case_AMI_comparison_M8_SNR15.png', 'Resolution', 200);

fig2 = figure('Name','PS Case Study Gain Comparison','Color','w');
bar(Gain_val);
grid on; box on;
set(gca, 'XTick', 1:numel(methods), 'XTickLabel', methods, 'XTickLabelRotation', 25);
ylabel('Validated gain over Uniform PAM (bits/symbol)');
title('M=8, SNR=15 dB: Gain comparison');
exportgraphics(fig2, 'PS_case_gain_comparison_M8_SNR15.png', 'Resolution', 200);

fprintf('[Plot] Saved PS_case_AMI_comparison_M8_SNR15.png\n');
fprintf('[Plot] Saved PS_case_gain_comparison_M8_SNR15.png\n');
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
    fprintf('    AMI fast    = %.8f bits/symbol\n', ami_fast);
    fprintf('    AMI val     = %.8f bits/symbol\n', ami_val);
    fprintf('    gain val    = %.8f bits/symbol\n\n', gain_val);
end


function print_x_p(x, p, P_avg)
    fprintf('   x = ['); fprintf(' %.6f', x(:)); fprintf(' ]\n');
    fprintf('   p = ['); fprintf(' %.6e', p(:)); fprintf(' ]\n');
    fprintf('   sum(p)=%.12f | dot(p,x)=%.12f | P_avg=%.12f | min(p)=%.3e\n\n', ...
        sum(p), dot(p(:),x(:)), P_avg, min(p));
end


function [x_best_final, info] = refine_gs_top_candidates(cfg, results_gs, px, params, ghN_h, y_grid, ami_pam_val, topK)
    info = struct();
    info.usedFmincon = false;
    info.candidates = [];

    allMI = [results_gs.bestMI];
    [~, idxSort] = sort(allMI, 'descend');
    topK = min(topK, numel(idxSort));

    % Start with best SA candidate by validated AMI.
    bestVal = -Inf;
    x_best_final = results_gs(idxSort(1)).x_best(:);

    for t = 1:topK
        idx = idxSort(t);
        x0 = AMI_functions.project_constellation_1D(results_gs(idx).x_best(:), cfg);
        val0 = AMI_functions.AMI_noCSI_validate(x0, px, params);
        if val0 > bestVal
            bestVal = val0;
            x_best_final = x0(:);
        end
    end

    if exist('fmincon', 'file') ~= 2
        fprintf('    fmincon unavailable. Keeping best SA candidate by validated AMI.\n');
        return;
    end

    fprintf('    Refining topK=%d GS candidates using fmincon...\n', topK);
    info.usedFmincon = true;

    M = cfg.M;
    minGap = cfg.SA.minGap;

    Aeq = ones(1,M) / M;
    beq = cfg.P_avg;
    lb = zeros(M,1);
    ub = [];

    A = zeros(M-1, M);
    b = -minGap * ones(M-1,1);
    for i = 1:M-1
        A(i,i)   = 1;
        A(i,i+1) = -1;
    end

    objective = @(x) -AMI_functions.AMI_noCSI_fast_grid(x(:).', px(:).', params, ghN_h, y_grid);

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
        x0 = AMI_functions.project_constellation_1D(results_gs(idx).x_best(:), cfg);

        try
            [x_ref, negMI, exitflag] = fmincon(objective, x0(:), A, b, Aeq, beq, lb, ub, [], options);
            x_ref = AMI_functions.project_constellation_1D(x_ref(:), cfg);
            amiFast = -negMI;
            amiVal = AMI_functions.AMI_noCSI_validate(x_ref, px, params);
        catch ME
            fprintf('        Refinement failed for candidate %d: %s\n', t, ME.message);
            x_ref = x0(:);
            exitflag = NaN;
            amiFast = AMI_functions.AMI_noCSI_fast_grid(x_ref(:).', px(:).', params, ghN_h, y_grid);
            amiVal = AMI_functions.AMI_noCSI_validate(x_ref, px, params);
        end

        cand(t).sourceStart = results_gs(idx).startIndex;
        cand(t).x = x_ref(:);
        cand(t).amiFast = amiFast;
        cand(t).amiVal = amiVal;
        cand(t).gainVal = amiVal - ami_pam_val;
        cand(t).exitflag = exitflag;

        fprintf('        [GS refine %d/%d] start #%d | AMI val=%.8f | gain=%.8f | exitflag=%g\n', ...
            t, topK, cand(t).sourceStart, cand(t).amiVal, cand(t).gainVal, cand(t).exitflag);

        if amiVal > bestVal
            bestVal = amiVal;
            x_best_final = x_ref(:);
        end
    end

    info.candidates = cand;
end
