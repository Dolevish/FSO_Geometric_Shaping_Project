function run_Joint_fmincon_refine_M8_SNR15()
%RUN_JOINT_FMINCON_REFINE_M8_SNR15
% Joint fmincon refinement for the best Alternating GS/PS solution.
%
% This script starts from the current best alternating solution for:
%   M=8, SNR=15 dB, sigma_X_sq=0.1, P_avg=1, minGap=0.005
%
% It then optimizes x and p together using fmincon:
%   maximize AMI_noCSI_fast_grid(x,p)
% subject to:
%   x_i >= 0
%   x_{i+1}-x_i >= minGap
%   p_i >= p_min
%   sum(p)=1
%   p'*x=P_avg
%
% Finally it validates the refined solution using AMI_noCSI_validate.

    clear; clc; close all;
    addpath(pwd);

    fprintf('============================================================\n');
    fprintf('   Joint fmincon refinement from Alternating GS/PS solution\n');
    fprintf('============================================================\n\n');

    %% ------------------------------------------------------------
    % 1. Parameters
    % -------------------------------------------------------------
    M            = 8;
    P_avg        = 1;
    SNR_dB       = 15;
    sigma_X_sq   = 0.1;
    R            = 1;
    minGap       = 0.005;
    p_min        = 1e-8;
    ghN_h        = 40;

    params = struct();
    params.M          = M;
    params.P_avg      = P_avg;
    params.SNR_dB     = SNR_dB;
    params.sigma_X_sq = sigma_X_sq;
    params.sigma_n_sq = P_avg^2 / (10^(SNR_dB/10));
    params.R          = R;

    sig_t_sq    = log(1 + sigma_X_sq);
    params.sig_t = sqrt(sig_t_sq);
    params.mu_t  = -0.5 * sig_t_sq;

    %% ------------------------------------------------------------
    % 2. Baseline and starting point
    % -------------------------------------------------------------
    x_pam = define_constellation(M, P_avg, true, "mean");
    p_uni = ones(M,1)/M;
    ami_pam_val = AMI_functions.AMI_noCSI_validate(x_pam, p_uni, params);

    % Current best alternating result from the 10-round run.
    % These values are used as a robust fallback even if the .mat file is not
    % available in the current MATLAB folder.
    x0 = [0.000000 0.154946 0.327682 0.476015 0.483059 0.885586 1.813304 4.373431].';
    p0 = [4.024365e-01 1.000000e-08 1.000000e-08 1.582905e-02 1.547552e-01 1.419546e-01 1.775456e-01 1.074790e-01].';

    % Make the fallback exactly feasible after rounding.
    p0 = project_probabilities_power(p0, x0, P_avg, p_min);

    % Try to load a more precise saved solution if it exists.
    if exist('result_Alternating_GS_PS_M8_SNR15.mat', 'file') == 2
        try
            S = load('result_Alternating_GS_PS_M8_SNR15.mat');
            [x_loaded, p_loaded, ok] = extract_best_xp_from_loaded_result(S, M);
            if ok
                x0 = x_loaded(:);
                p0 = p_loaded(:);
                p0 = project_probabilities_power(p0, x0, P_avg, p_min);
                fprintf('[Load] Loaded x0,p0 from result_Alternating_GS_PS_M8_SNR15.mat\n');
            else
                fprintf('[Load] MAT file found, but best x/p fields were not recognized. Using hardcoded latest solution.\n');
            end
        catch ME
            fprintf('[Load] Could not load MAT result safely: %s\n', ME.message);
            fprintf('[Load] Using hardcoded latest solution.\n');
        end
    else
        fprintf('[Load] result_Alternating_GS_PS_M8_SNR15.mat not found. Using hardcoded latest solution.\n');
    end

    % Build a wide grid for initial reporting and fmincon objective.
    xMaxUB = 8;
    y_grid = AMI_functions.build_noCSI_y_grid(params, xMaxUB);

    ami0_fast = AMI_functions.AMI_noCSI_fast_grid(x0(:).', p0(:).', params, ghN_h, y_grid);
    ami0_val  = AMI_functions.AMI_noCSI_validate(x0(:).', p0(:).', params);

    fprintf('\n[1] Case parameters\n');
    fprintf('    M              = %d\n', M);
    fprintf('    SNR            = %.1f dB\n', SNR_dB);
    fprintf('    P_avg          = %.6f\n', P_avg);
    fprintf('    sigma_X_sq     = %.6f\n', sigma_X_sq);
    fprintf('    sigma_n_sq     = %.6e\n', params.sigma_n_sq);
    fprintf('    minGap         = %.6f\n', minGap);
    fprintf('    p_min          = %.2e\n', p_min);
    fprintf('    xMaxUB         = %.3f\n', xMaxUB);
    fprintf('    y_grid range   = [%.4f, %.4f], Ny=%d\n', y_grid(1), y_grid(end), numel(y_grid));

    fprintf('\n[2] Starting point: best Alternating GS/PS\n');
    print_solution('Initial alternating solution', x0, p0, P_avg, minGap, ami0_fast, ami0_val, ami_pam_val);

    %% ------------------------------------------------------------
    % 3. Joint fmincon refinement
    % -------------------------------------------------------------
    opts = struct();
    opts.minGap       = minGap;
    opts.p_min        = p_min;
    opts.ghN_h        = ghN_h;
    opts.xMaxUB       = xMaxUB;
    opts.y_grid       = y_grid;

    % A few starts around the current best solution. Start #1 is exactly x0,p0.
    opts.nStarts      = 5;
    opts.xPerturbStd  = 0.02;
    opts.pPerturbStd  = 0.01;

    opts.maxIter      = 600;
    opts.maxFunEvals  = 30000;
    opts.display      = 'iter';
    opts.seed         = 20260527;

    % Usually keep this false. Setting true may help only if Parallel Toolbox
    % is available and the objective is expensive enough.
    opts.useParallelFiniteDiff = false;

    out_joint = joint_fmincon_refine_gs_ps(x0, p0, params, opts);

    %% ------------------------------------------------------------
    % 4. Final comparison
    % -------------------------------------------------------------
    xj = out_joint.x_best(:);
    pj = out_joint.p_best(:);

    ami_joint_fast = out_joint.ami_fast;
    ami_joint_val  = out_joint.ami_val;

    fprintf('\n============================================================\n');
    fprintf(' Final comparison\n');
    fprintf('============================================================\n');
    fprintf('--------------------------------------------------------------------------\n');
    fprintf(' Method                       | AMI fast   | AMI validated | Gain val\n');
    fprintf('--------------------------------------------------------------------------\n');
    fprintf(' Uniform PAM                  | %10.8f | %13.8f | %10.8f\n', ...
        NaN, ami_pam_val, 0);
    fprintf(' Alternating GS/PS start      | %10.8f | %13.8f | %10.8f\n', ...
        ami0_fast, ami0_val, ami0_val - ami_pam_val);
    fprintf(' Joint fmincon refined        | %10.8f | %13.8f | %10.8f\n', ...
        ami_joint_fast, ami_joint_val, ami_joint_val - ami_pam_val);
    fprintf('--------------------------------------------------------------------------\n');
    fprintf(' Improvement over start, validated AMI = %.10f bits/symbol\n', ami_joint_val - ami0_val);

    print_solution('Joint fmincon refined solution', xj, pj, P_avg, minGap, ami_joint_fast, ami_joint_val, ami_pam_val);

    %% ------------------------------------------------------------
    % 5. Save and plot
    % -------------------------------------------------------------
    result = struct();
    result.params        = params;
    result.opts          = opts;
    result.x_pam         = x_pam;
    result.p_uni         = p_uni;
    result.ami_pam_val   = ami_pam_val;
    result.x_start       = x0;
    result.p_start       = p0;
    result.ami_start_fast= ami0_fast;
    result.ami_start_val = ami0_val;
    result.out_joint     = out_joint;
    result.x_joint       = xj;
    result.p_joint       = pj;
    result.ami_joint_fast= ami_joint_fast;
    result.ami_joint_val = ami_joint_val;
    result.gain_joint_val= ami_joint_val - ami_pam_val;

    save('result_Joint_fmincon_refine_M8_SNR15.mat', 'result');
    fprintf('\n[Save] Saved result_Joint_fmincon_refine_M8_SNR15.mat\n');

    fig = figure('Color','w', 'Name','Joint fmincon GS/PS refinement');
    tiledlayout(2,1);

    nexttile;
    stem(1:M, x0, 'o', 'LineWidth', 1.4); hold on;
    stem(1:M, xj, 'x', 'LineWidth', 1.4);
    grid on; box on;
    xlabel('Symbol index'); ylabel('x_i');
    title('Constellation before/after joint fmincon');
    legend('Alternating start','Joint refined','Location','best');

    nexttile;
    stem(1:M, p0, 'o', 'LineWidth', 1.4); hold on;
    stem(1:M, pj, 'x', 'LineWidth', 1.4);
    grid on; box on;
    xlabel('Symbol index'); ylabel('p_i');
    title('Probabilities before/after joint fmincon');
    legend('Alternating start','Joint refined','Location','best');

    exportgraphics(fig, 'Joint_fmincon_refine_M8_SNR15_x_p.png', 'Resolution', 300);
    fprintf('[Plot] Saved Joint_fmincon_refine_M8_SNR15_x_p.png\n');

    fprintf('\nDone.\n');
end

% =========================================================================
% Helpers
% =========================================================================
function print_solution(name, x, p, P_avg, minGap, ami_fast, ami_val, ami_pam_val)
    fprintf('\n    %s\n', name);
    fprintf('    x = ['); fprintf(' %.6f', x); fprintf(' ]\n');
    fprintf('    p = ['); fprintf(' %.6e', p); fprintf(' ]\n');
    fprintf('    sum(p)      = %.12f\n', sum(p));
    fprintf('    dot(p,x)    = %.12f\n', dot(p,x));
    fprintf('    target power= %.12f\n', P_avg);
    fprintf('    min(p)      = %.3e\n', min(p));
    fprintf('    min diff(x) = %.3e\n', min(diff(x)));
    fprintf('    minGap      = %.3e\n', minGap);
    fprintf('    AMI fast    = %.8f bits/symbol\n', ami_fast);
    fprintf('    AMI val     = %.8f bits/symbol\n', ami_val);
    fprintf('    gain val    = %.8f bits/symbol\n', ami_val - ami_pam_val);
end

function [x, p, ok] = extract_best_xp_from_loaded_result(S, M)
    ok = false;
    x = [];
    p = [];

    % Most likely format from our generated script.
    if isfield(S, 'result')
        R = S.result;
        candidates = {};

        if isfield(R, 'x_best') && isfield(R, 'p_best')
            candidates{end+1} = {R.x_best, R.p_best}; 
        end
        if isfield(R, 'best') && isstruct(R.best) && isfield(R.best,'x') && isfield(R.best,'p')
            candidates{end+1} = {R.best.x, R.best.p}; 
        end
        if isfield(R, 'bestResult') && isstruct(R.bestResult) && isfield(R.bestResult,'x') && isfield(R.bestResult,'p')
            candidates{end+1} = {R.bestResult.x, R.bestResult.p}; 
        end
        if isfield(R, 'history')
            H = R.history;
            try
                if isstruct(H)
                    vals = [H.ami_val];
                    [~, idx] = max(vals);
                    candidates{end+1} = {H(idx).x, H(idx).p}; 
                end
            catch
            end
        end

        for k = 1:numel(candidates)
            xt = candidates{k}{1};
            pt = candidates{k}{2};
            if numel(xt) == M && numel(pt) == M
                x = xt(:);
                p = pt(:);
                ok = true;
                return;
            end
        end
    end

    % Fallback: search top-level fields.
    if isfield(S, 'x_best') && isfield(S, 'p_best') && numel(S.x_best)==M && numel(S.p_best)==M
        x = S.x_best(:);
        p = S.p_best(:);
        ok = true;
    end
end
