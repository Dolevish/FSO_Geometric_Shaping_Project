%% run_GS_M8_SNR_sweep_minGap0.m
% Comprehensive GS simulation for M=8 with no minimum-spacing constraint.
%
% Case:
%   M = 8
%   SNR = 0:5:30 dB
%   minGap = 0
%   sigma_R^2 = 0.1 by default
%
% Goal:
%   For each SNR value:
%       1. Build uniform positive PAM baseline.
%       2. Compute validated PAM AMI.
%       3. Run multi-start simulated annealing for geometric shaping.
%       4. Run local constrained refinement with fmincon on the top SA candidates.
%       5. Select the best final candidate by validated AMI.
%       6. Print the final constellation, AMI, and gain.
%   Finally:
%       Plot validated AMI gain versus SNR.
%
% Notes:
%   - minGap = 0 intentionally allows repeated/merged constellation levels.
%   - This can produce degenerate GS solutions, which are useful as a bridge
%     toward later comparison with probabilistic shaping (PS).
%   - Requires the project files:
%       AMI_functions.m
%       calculate_Py_given_x.m
%       define_constellation.m
%       Enforce_Power_Constraint.m
%       simulated_annealing.m
%       sa_multistart.m
%   - fmincon is optional. If Optimization Toolbox is unavailable, the script
%     uses the best SA-only result.

clear; clc; close all;
addpath(pwd);

fprintf('============================================================\n');
fprintf('   GS SNR Sweep: M=8, minGap=0, SA + Local Refinement\n');
fprintf('============================================================\n\n');

%% 0. Sweep definition

SNR_dB_list = 0:5:30;
numSNR = numel(SNR_dB_list);

fprintf('[0] SNR sweep values\n');
fprintf('    SNR_dB_list = ['); fprintf(' %.1f', SNR_dB_list); fprintf(' ]\n\n');

%% 1. Fixed simulation settings

baseParams = struct();
baseParams.R          = 1;      % Detector responsivity
baseParams.M          = 8;      % Constellation size
baseParams.P_avg      = 1;      % Average optical power
baseParams.sigma_X_sq = 0.1;    % sigma_R^2 turbulence parameter; change if needed

imdd_mode       = true;
powerConstraint = "mean";

minGap = 0;
ghN_h  = 40;

% SA settings. Adjust these if runtime is too long/short.
SA_template = struct();
SA_template.imdd_mode       = true;
SA_template.enforce_sort    = true;
SA_template.enforce_power   = true;
SA_template.powerConstraint = "mean";

SA_template.minGap       = minGap;
SA_template.projectIters = 3;

SA_template.nStarts       = 12;
SA_template.maxIter       = 8000;
SA_template.itersPerTemp  = 100;

SA_template.T0            = 0.4;
SA_template.Tf            = 1e-3;
SA_template.nBlocks       = ceil(SA_template.maxIter / SA_template.itersPerTemp);
SA_template.coolingRate   = exp(log(SA_template.Tf / SA_template.T0) / SA_template.nBlocks);

SA_template.baseStd0      = 0.20;
SA_template.baseStdMin    = 1e-4;
SA_template.baseStdMax    = 0.80;
SA_template.baseStdGrow   = 1.10;
SA_template.baseStdShrink = 0.80;

SA_template.targetAccLo   = 0.20;
SA_template.targetAccHi   = 0.60;

SA_template.logEvery      = 0;

% Use [] for random seed each full script run.
% For reproducible SNR-specific runs, this script derives a seed from this base.
SA_template.seedInit      = 12345;

SA_template.useParallel       = true;
SA_template.numWorkers        = [];
SA_template.closePoolWhenDone = false;

% fmincon local refinement settings
refine = struct();
refine.enabled = true;
refine.topK    = 5;
refine.maxIterations = 300;
refine.maxFunctionEvaluations = 5000;
refine.algorithm = 'sqp';

fprintf('[1] Fixed settings\n');
fprintf('    M              = %d\n', baseParams.M);
fprintf('    P_avg          = %.4f\n', baseParams.P_avg);
fprintf('    sigma_X_sq     = %.4f\n', baseParams.sigma_X_sq);
fprintf('    minGap         = %.4f\n', minGap);
fprintf('    ghN_h          = %d\n', ghN_h);
fprintf('    nStarts        = %d\n', SA_template.nStarts);
fprintf('    maxIter        = %d\n', SA_template.maxIter);
fprintf('    refinement topK= %d\n\n', refine.topK);

%% 2. Check dependencies

requiredFiles = { ...
    'AMI_functions.m', ...
    'calculate_Py_given_x.m', ...
    'define_constellation.m', ...
    'Enforce_Power_Constraint.m', ...
    'simulated_annealing.m', ...
    'sa_multistart.m'};

for k = 1:numel(requiredFiles)
    if ~exist(requiredFiles{k}, 'file')
        error('Required file not found in path: %s', requiredFiles{k});
    end
end

hasFmincon = exist('fmincon', 'file') == 2;
if refine.enabled && ~hasFmincon
    warning('fmincon not found. Local refinement will be skipped.');
    refine.enabled = false;
end

%% 3. Preallocate results

sweep = repmat(struct( ...
    'SNR_dB', NaN, ...
    'params', [], ...
    'x_pam', [], ...
    'ami_pam_fast', NaN, ...
    'ami_pam_val', NaN, ...
    'x_sa_best', [], ...
    'ami_sa_fast', NaN, ...
    'ami_sa_val', NaN, ...
    'gain_sa_val', NaN, ...
    'x_refined_best', [], ...
    'ami_refined_fast', NaN, ...
    'ami_refined_val', NaN, ...
    'gain_refined_val', NaN, ...
    'x_final', [], ...
    'ami_final_fast', NaN, ...
    'ami_final_val', NaN, ...
    'gain_final_fast', NaN, ...
    'gain_final_val', NaN, ...
    'selectedSource', '', ...
    'out_sa', [], ...
    'results_sa', [], ...
    'refineResults', []), numSNR, 1);

%% 4. Main SNR loop

for sIdx = 1:numSNR
    SNR_dB = SNR_dB_list(sIdx);

    fprintf('\n============================================================\n');
    fprintf('   SNR case %d/%d: SNR = %.1f dB\n', sIdx, numSNR, SNR_dB);
    fprintf('============================================================\n\n');

    %% A. Build channel parameters for this SNR

    params = baseParams;
    params.SNR_dB = SNR_dB;

    % Paper convention: SNR = P_avg^2 / sigma_n^2
    params.sigma_n_sq = params.P_avg^2 / 10^(params.SNR_dB/10);

    sig_t_sq     = log(1 + params.sigma_X_sq);
    params.sig_t = sqrt(sig_t_sq);
    params.mu_t  = -0.5 * sig_t_sq;

    fprintf('[A] Parameters\n');
    fprintf('    SNR           = %.1f dB\n', params.SNR_dB);
    fprintf('    sigma_n_sq    = %.6e\n', params.sigma_n_sq);
    fprintf('    sigma_X_sq    = %.6f\n', params.sigma_X_sq);
    fprintf('    mu_t          = %.6f\n', params.mu_t);
    fprintf('    sig_t         = %.6f\n\n', params.sig_t);

    %% B. Baseline uniform positive PAM

    x_pam = define_constellation(params.M, params.P_avg, imdd_mode, powerConstraint);
    px = ones(params.M, 1) / params.M;

    fprintf('[B] Uniform PAM baseline\n');
    fprintf('    x_pam = ['); fprintf(' %.6f', x_pam); fprintf(' ]\n');
    fprintf('    mean(x_pam) = %.12f\n\n', mean(x_pam));

    %% C. Build y-grid and configure AMI evaluator

    % Conservative y bound. At minGap=0, optimized solutions can push the last
    % point above the uniform PAM maximum, so use a wide bound.
    x_max_bound = params.M * params.P_avg * 1.2; % Gives 9.6, a very safe margin
    y_grid = AMI_functions.build_noCSI_y_grid(params, x_max_bound);

    fprintf('[C] Fast AMI grid\n');
    fprintf('    y range = [%.4f, %.4f]\n', y_grid(1), y_grid(end));
    fprintf('    Ny      = %d\n\n', numel(y_grid));

    cfg = struct();
    cfg.params = params;
    cfg.M      = params.M;
    cfg.P_avg  = params.P_avg;
    cfg.SA     = SA_template;

    % Use a different reproducible seed per SNR so cases do not repeat exactly.
    if ~isempty(SA_template.seedInit)
        cfg.SA.seedInit = SA_template.seedInit + 1000*sIdx;
    end

    cfg.AMI_Evaluator = @(x) AMI_functions.AMI_noCSI_fast_grid( ...
        x(:).', px, params, ghN_h, y_grid);

    %% D. Baseline AMI

    fprintf('[D] Computing PAM AMI...\n');
    ami_pam_fast = AMI_functions.AMI_noCSI_fast_grid(x_pam, px, params, ghN_h, y_grid);
    ami_pam_val  = AMI_functions.AMI_noCSI_validate(x_pam, px, params);

    fprintf('    PAM AMI fast      = %.8f bits/symbol\n', ami_pam_fast);
    fprintf('    PAM AMI validated = %.8f bits/symbol\n\n', ami_pam_val);

    %% E. Run multi-start SA

    fprintf('[E] Running SA multistart for SNR=%.1f dB...\n\n', SNR_dB);
    taskLabel = sprintf('[M8-SNR%02.0f-gap0]', SNR_dB);
    [out_sa, results_sa] = sa_multistart(cfg, x_pam(:), taskLabel);

    x_sa_best = AMI_functions.project_constellation_1D(out_sa.bestX(:), cfg);
    x_sa_best = x_sa_best(:).';

    ami_sa_fast = AMI_functions.AMI_noCSI_fast_grid(x_sa_best, px, params, ghN_h, y_grid);
    ami_sa_val  = AMI_functions.AMI_noCSI_validate(x_sa_best, px, params);
    gain_sa_val = ami_sa_val - ami_pam_val;

    fprintf('\n[E] Best SA-only result\n');
    fprintf('    best start       = %d\n', out_sa.bestStart);
    fprintf('    x_sa_best        = ['); fprintf(' %.6f', x_sa_best); fprintf(' ]\n');
    fprintf('    SA AMI fast      = %.8f bits/symbol\n', ami_sa_fast);
    fprintf('    SA AMI validated = %.8f bits/symbol\n', ami_sa_val);
    fprintf('    SA gain validated= %.8f bits/symbol\n\n', gain_sa_val);

    %% F. Local constrained refinement on top SA candidates

    refineResults = [];
    x_refined_best = [];
    ami_refined_fast = -Inf;
    ami_refined_val  = -Inf;
    gain_refined_val = -Inf;

    if refine.enabled
        fprintf('[F] Running local constrained refinement on top SA candidates...\n');

        [~, idxSort] = sort([results_sa.bestMI], 'descend');
        topK = min(refine.topK, numel(idxSort));
        idxTop = idxSort(1:topK);

        fprintf('    Refining topK = %d SA candidates using fmincon...\n\n', topK);

        refineResults = repmat(struct( ...
            'rank', NaN, ...
            'sourceStart', NaN, ...
            'x0', [], ...
            'x_refined', [], ...
            'ami_fast', NaN, ...
            'ami_val', NaN, ...
            'gain_val', NaN, ...
            'exitflag', NaN, ...
            'output', []), topK, 1);

        for r = 1:topK
            srcIdx = idxTop(r);
            x0_refine = results_sa(srcIdx).x_best(:);
            x0_refine = AMI_functions.project_constellation_1D(x0_refine, cfg);
            x0_refine = x0_refine(:);

            fprintf('    [Refine %d/%d] From SA start #%d | SA fast AMI = %.8f\n', ...
                r, topK, srcIdx, results_sa(srcIdx).bestMI);

            [x_ref, ami_ref_fast, exitflag, output] = refine_with_fmincon( ...
                x0_refine, cfg, px, params, ghN_h, y_grid, refine);

            x_ref = AMI_functions.project_constellation_1D(x_ref(:), cfg);
            x_ref = x_ref(:).';

            ami_ref_fast = AMI_functions.AMI_noCSI_fast_grid(x_ref, px, params, ghN_h, y_grid);
            ami_ref_val  = AMI_functions.AMI_noCSI_validate(x_ref, px, params);
            gain_ref_val = ami_ref_val - ami_pam_val;

            refineResults(r).rank        = r;
            refineResults(r).sourceStart = srcIdx;
            refineResults(r).x0          = x0_refine(:).';
            refineResults(r).x_refined   = x_ref;
            refineResults(r).ami_fast    = ami_ref_fast;
            refineResults(r).ami_val     = ami_ref_val;
            refineResults(r).gain_val    = gain_ref_val;
            refineResults(r).exitflag    = exitflag;
            refineResults(r).output      = output;

            fprintf('        Refined x          = ['); fprintf(' %.6f', x_ref); fprintf(' ]\n');
            fprintf('        Refined AMI fast   = %.8f bits/symbol\n', ami_ref_fast);
            fprintf('        Refined AMI val    = %.8f bits/symbol\n', ami_ref_val);
            fprintf('        Refined gain val   = %.8f bits/symbol\n', gain_ref_val);
            fprintf('        exitflag           = %d\n\n', exitflag);
        end

        [ami_refined_val, bestR] = max([refineResults.ami_val]);
        x_refined_best  = refineResults(bestR).x_refined;
        ami_refined_fast = refineResults(bestR).ami_fast;
        gain_refined_val = refineResults(bestR).gain_val;

        fprintf('[F] Best refined candidate selected by validated AMI\n');
        fprintf('    source SA start  = %d\n', refineResults(bestR).sourceStart);
        fprintf('    best refined x   = ['); fprintf(' %.6f', x_refined_best); fprintf(' ]\n');
        fprintf('    best AMI fast    = %.8f bits/symbol\n', ami_refined_fast);
        fprintf('    best AMI val     = %.8f bits/symbol\n\n', ami_refined_val);
    else
        fprintf('[F] Local refinement skipped.\n\n');
    end

    %% G. Select final candidate by validated AMI

    if refine.enabled && isfinite(ami_refined_val) && ami_refined_val >= ami_sa_val
        x_final        = x_refined_best;
        ami_final_fast = ami_refined_fast;
        ami_final_val  = ami_refined_val;
        selectedSource = 'Hybrid refined';
    else
        x_final        = x_sa_best;
        ami_final_fast = ami_sa_fast;
        ami_final_val  = ami_sa_val;
        selectedSource = 'SA only';
    end

    gain_final_fast = ami_final_fast - ami_pam_fast;
    gain_final_val  = ami_final_val  - ami_pam_val;

    fprintf('[G] Final validation and gain report for SNR=%.1f dB\n', SNR_dB);
    fprintf('    selected source   = %s\n', selectedSource);
    fprintf('    Final GS x        = ['); fprintf(' %.6f', x_final); fprintf(' ]\n');
    fprintf('    mean(x_final)     = %.12f\n', mean(x_final));
    fprintf('    min(x_final)      = %.12f\n', min(x_final));
    fprintf('    min diff(x_final) = %.12f\n\n', min(diff(x_final)));

    fprintf('    PAM AMI fast      = %.8f bits/symbol\n', ami_pam_fast);
    fprintf('    PAM AMI validated = %.8f bits/symbol\n', ami_pam_val);
    fprintf('    GS AMI fast       = %.8f bits/symbol\n', ami_final_fast);
    fprintf('    GS AMI validated  = %.8f bits/symbol\n', ami_final_val);
    fprintf('    AMI gain fast     = %.8f bits/symbol\n', gain_final_fast);
    fprintf('    AMI gain validated= %.8f bits/symbol\n\n', gain_final_val);

    %% H. Store results

    sweep(sIdx).SNR_dB            = SNR_dB;
    sweep(sIdx).params            = params;
    sweep(sIdx).x_pam             = x_pam;
    sweep(sIdx).ami_pam_fast      = ami_pam_fast;
    sweep(sIdx).ami_pam_val       = ami_pam_val;
    sweep(sIdx).x_sa_best         = x_sa_best;
    sweep(sIdx).ami_sa_fast       = ami_sa_fast;
    sweep(sIdx).ami_sa_val        = ami_sa_val;
    sweep(sIdx).gain_sa_val       = gain_sa_val;
    sweep(sIdx).x_refined_best    = x_refined_best;
    sweep(sIdx).ami_refined_fast  = ami_refined_fast;
    sweep(sIdx).ami_refined_val   = ami_refined_val;
    sweep(sIdx).gain_refined_val  = gain_refined_val;
    sweep(sIdx).x_final           = x_final;
    sweep(sIdx).ami_final_fast    = ami_final_fast;
    sweep(sIdx).ami_final_val     = ami_final_val;
    sweep(sIdx).gain_final_fast   = gain_final_fast;
    sweep(sIdx).gain_final_val    = gain_final_val;
    sweep(sIdx).selectedSource    = selectedSource;
    sweep(sIdx).out_sa            = out_sa;
    sweep(sIdx).results_sa        = results_sa;
    sweep(sIdx).refineResults     = refineResults;

    % Save checkpoint after every SNR in case the full run is interrupted.
    save('result_GS_M8_SNR_sweep_minGap0_checkpoint.mat', ...
        'sweep', 'SNR_dB_list', 'baseParams', 'SA_template', 'refine', 'ghN_h', 'minGap');
end

%% 5. Summary table

fprintf('\n============================================================\n');
fprintf('   SNR sweep summary: M=8, minGap=0\n');
fprintf('============================================================\n');
fprintf('---------------------------------------------------------------------------------------------\n');
fprintf(' idx | SNR[dB] | PAM val AMI | GS val AMI | Gain val | min diff | source\n');
fprintf('---------------------------------------------------------------------------------------------\n');

for sIdx = 1:numSNR
    fprintf(' %3d | %7.1f | %11.8f | %10.8f | %8.8f | %8.3e | %s\n', ...
        sIdx, ...
        sweep(sIdx).SNR_dB, ...
        sweep(sIdx).ami_pam_val, ...
        sweep(sIdx).ami_final_val, ...
        sweep(sIdx).gain_final_val, ...
        min(diff(sweep(sIdx).x_final)), ...
        sweep(sIdx).selectedSource);
end
fprintf('---------------------------------------------------------------------------------------------\n\n');

fprintf('Final constellations by SNR:\n');
for sIdx = 1:numSNR
    fprintf('  SNR = %4.1f dB | x = [', sweep(sIdx).SNR_dB);
    fprintf(' %.6f', sweep(sIdx).x_final);
    fprintf(' ] | AMI = %.8f | Gain = %.8f\n', ...
        sweep(sIdx).ami_final_val, sweep(sIdx).gain_final_val);
end
fprintf('\n');

%% 6. Save final results

result = struct();
result.sweep       = sweep;
result.SNR_dB_list = SNR_dB_list;
result.baseParams  = baseParams;
result.SA_template = SA_template;
result.refine      = refine;
result.ghN_h       = ghN_h;
result.minGap      = minGap;

save('result_GS_M8_SNR_sweep_minGap0.mat', 'result');
fprintf('[Save] Saved full results to result_GS_M8_SNR_sweep_minGap0.mat\n');

%% 7. Plot gain versus SNR

snr_vec      = [sweep.SNR_dB];
gain_val_vec = [sweep.gain_final_val];
gain_sa_vec  = [sweep.gain_sa_val];
ami_pam_vec  = [sweep.ami_pam_val];
ami_gs_vec   = [sweep.ami_final_val];

figure('Name', 'GS AMI Gain vs SNR, M=8, minGap=0');
plot(snr_vec, gain_val_vec, '-o', 'LineWidth', 1.8, 'MarkerSize', 7);
grid on; box on;
xlabel('SNR [dB]');
ylabel('Validated AMI gain [bits/symbol]');
title(sprintf('GS AMI Gain vs SNR, M=8, minGap=0, \\sigma_R^2=%.2f', baseParams.sigma_X_sq));

saveas(gcf, 'GS_gain_vs_SNR_M8_minGap0.png');
fprintf('[Plot] Saved gain plot to GS_gain_vs_SNR_M8_minGap0.png\n');

%% 8. Optional AMI comparison plot

figure('Name', 'PAM AMI vs GS AMI, M=8, minGap=0');
plot(snr_vec, ami_pam_vec, '-o', 'LineWidth', 1.8, 'MarkerSize', 7); hold on;
plot(snr_vec, ami_gs_vec, '-x', 'LineWidth', 1.8, 'MarkerSize', 8);
grid on; box on;
xlabel('SNR [dB]');
ylabel('Validated AMI [bits/symbol]');
title(sprintf('Validated AMI vs SNR, M=8, minGap=0, \\sigma_R^2=%.2f', baseParams.sigma_X_sq));
legend('Uniform PAM', 'GS minGap=0', 'Location', 'best');

saveas(gcf, 'AMI_vs_SNR_M8_minGap0.png');
fprintf('[Plot] Saved AMI comparison plot to AMI_vs_SNR_M8_minGap0.png\n');

fprintf('\nDone.\n');

%% ========================================================================
% Local helper function
% ========================================================================

function [x_refined, ami_fast, exitflag, output] = refine_with_fmincon( ...
    x0, cfg, px, params, ghN_h, y_grid, refine)
%REFINE_WITH_FMINCON Local constrained maximization of fast AMI.
%   fmincon minimizes, so the objective is negative AMI.
%
% Constraints:
%   x_i >= 0
%   mean(x) = P_avg
%   x_{i+1} - x_i >= minGap
% For minGap = 0, this simply enforces monotonic nondecreasing order.

    M = params.M;
    minGap = cfg.SA.minGap;

    objective = @(x) -AMI_functions.AMI_noCSI_fast_grid( ...
        x(:).', px, params, ghN_h, y_grid);

    % Equality: mean(x) = P_avg
    Aeq = ones(1, M) / M;
    beq = params.P_avg;

    % Inequality: x_i - x_{i+1} <= -minGap
    A = zeros(M-1, M);
    b = -minGap * ones(M-1, 1);
    for i = 1:(M-1)
        A(i,i)   = 1;
        A(i,i+1) = -1;
    end

    % Non-negativity
    lb = zeros(M, 1);
    ub = [];

    x0 = AMI_functions.project_constellation_1D(x0(:), cfg);
    x0 = x0(:);

    options = optimoptions('fmincon', ...
        'Display', 'off', ...
        'Algorithm', refine.algorithm, ...
        'MaxIterations', refine.maxIterations, ...
        'MaxFunctionEvaluations', refine.maxFunctionEvaluations, ...
        'OptimalityTolerance', 1e-8, ...
        'StepTolerance', 1e-10, ...
        'ConstraintTolerance', 1e-10);

    [x_refined, negAMI, exitflag, output] = fmincon( ...
        objective, x0, A, b, Aeq, beq, lb, ub, [], options);

    x_refined = AMI_functions.project_constellation_1D(x_refined(:), cfg);
    ami_fast  = -negAMI;
end
