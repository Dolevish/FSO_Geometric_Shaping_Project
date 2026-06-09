function validate_joint_solution_M8_SNR15_strict()
%VALIDATE_JOINT_SOLUTION_M8_SNR15_STRICT
% Strict validation for the previous best joint GS+PS solution.
%
% This script validates the solution found in the previous joint fmincon run:
%   x = [0, 0.179700, 0.411914, 0.521979, 1.079611, 1.920520, 3.841593, 7.666677]
%   p = [0.4215765, 1e-8, 1e-8, 0.2104800, 0.1336925, 0.1278013, 0.0825505, 0.0238992]
%
% It performs several independent checks:
%   1. Standard AMI_functions.AMI_noCSI_validate.
%   2. Strict adaptive-Pyx validation on wider and denser custom y-grids.
%   3. PDF mass checks: integral p(y|x_i) dy should be close to 1.
%   4. Fast GH/grid sensitivity check for ghN_h = 40, 60, 80.
%
% Required files:
%   AMI_functions.m
%   calculate_Py_given_x.m
%   define_constellation.m
%   project_probabilities_power.m

    clear; clc; close all;

    fprintf('============================================================\n');
    fprintf('   Strict validation: final joint GS/PS solution, M=8 SNR=15\n');
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

    %% ------------------------------------------------------------
    % 2. Previous best solution from joint fmincon
    % -------------------------------------------------------------
    x = [0.000000 0.179700 0.411914 0.521979 1.079611 1.920520 3.841593 7.666677].';
    p = [4.215765e-01 1.000000e-08 1.000000e-08 2.104800e-01 1.336925e-01 1.278013e-01 8.255049e-02 2.389924e-02].';

    % Repair p after decimal rounding while keeping x fixed.
    p = project_probabilities_power(p, x, P_avg, p_min);

    x_pam = define_constellation(M, P_avg, true, "mean");
    p_uni = ones(M,1)/M;
    ami_pam_val = AMI_functions.AMI_noCSI_validate(x_pam(:).', p_uni(:).', params);

    fprintf('[1] Case parameters\n');
    fprintf('    M              = %d\n', M);
    fprintf('    SNR            = %.1f dB\n', SNR_dB);
    fprintf('    P_avg          = %.6f\n', P_avg);
    fprintf('    sigma_X_sq     = %.6f\n', sigma_X_sq);
    fprintf('    sigma_n_sq     = %.6e\n', params.sigma_n_sq);
    fprintf('    mu_t           = %.6f\n', params.mu_t);
    fprintf('    sig_t          = %.6f\n', params.sig_t);

    fprintf('\n[2] Solution to validate\n');
    print_constraints(x, p, P_avg, minGap);
    fprintf('    x = ['); fprintf(' %.6f', x); fprintf(' ]\n');
    fprintf('    p = ['); fprintf(' %.6e', p); fprintf(' ]\n');

    %% ------------------------------------------------------------
    % 3. Standard validation
    % -------------------------------------------------------------
    fprintf('\n[3] Standard validation using AMI_functions.AMI_noCSI_validate\n');
    t0 = tic;
    ami_standard = AMI_functions.AMI_noCSI_validate(x(:).', p(:).', params);
    rt_standard = toc(t0);
    fprintf('    AMI standard validated = %.10f bits/symbol | runtime %.2f sec\n', ami_standard, rt_standard);
    fprintf('    Gain vs PAM            = %.10f bits/symbol\n', ami_standard - ami_pam_val);

    %% ------------------------------------------------------------
    % 4. Fast GH/grid sensitivity
    % -------------------------------------------------------------
    fprintf('\n[4] Fast evaluator sensitivity check\n');
    gh_list = [40 60 80];
    xMaxBounds = [8 10 12];
    fastTable = [];

    fprintf('    ghN | xMaxBound | Ny    | AMI fast\n');
    fprintf('    ----+-----------+-------+-------------\n');
    for a = 1:numel(gh_list)
        for b = 1:numel(xMaxBounds)
            ghN = gh_list(a);
            xb  = xMaxBounds(b);
            yg  = AMI_functions.build_noCSI_y_grid(params, xb);
            mi_fast = AMI_functions.AMI_noCSI_fast_grid(x(:).', p(:).', params, ghN, yg);
            fastTable = [fastTable; ghN, xb, numel(yg), mi_fast]; %#ok<AGROW>
            fprintf('    %3d | %9.2f | %5d | %.10f\n', ghN, xb, numel(yg), mi_fast);
        end
    end

    %% ------------------------------------------------------------
    % 5. Strict adaptive-Pyx validation on custom y-grids
    % -------------------------------------------------------------
    fprintf('\n[5] Strict adaptive-Pyx validation on custom y-grids\n');

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

    strictResults = repmat(struct('name','','Ny',NaN,'yMin',NaN,'yMax',NaN, ...
        'ami',NaN,'gain',NaN,'runtime',NaN,'massMin',NaN,'massMax',NaN, ...
        'massErrMax',NaN,'masses',[]), numel(configs), 1);

    fprintf('    config          | Ny     | y-range              | AMI strict  | gain       | max PDF mass err | runtime\n');
    fprintf('    ----------------+--------+----------------------+-------------+------------+------------------+---------\n');

    for k = 1:numel(configs)
        yg = build_strict_y_grid(params, x, configs(k));
        t1 = tic;
        [mi_strict, masses] = ami_validate_on_custom_grid(x, p, params, yg);
        rt = toc(t1);
        massErrMax = max(abs(masses - 1));

        strictResults(k).name = configs(k).name;
        strictResults(k).Ny = numel(yg);
        strictResults(k).yMin = yg(1);
        strictResults(k).yMax = yg(end);
        strictResults(k).ami = mi_strict;
        strictResults(k).gain = mi_strict - ami_pam_val;
        strictResults(k).runtime = rt;
        strictResults(k).massMin = min(masses);
        strictResults(k).massMax = max(masses);
        strictResults(k).massErrMax = massErrMax;
        strictResults(k).masses = masses;

        fprintf('    %-15s | %6d | [%8.3f,%8.3f] | %.10f | %.10f | %.3e        | %.1fs\n', ...
            configs(k).name, numel(yg), yg(1), yg(end), mi_strict, mi_strict - ami_pam_val, massErrMax, rt);
    end

    %% ------------------------------------------------------------
    % 6. Final decision summary
    % -------------------------------------------------------------
    strictAMIs = [strictResults.ami];
    ami_min = min([ami_standard, strictAMIs]);
    ami_max = max([ami_standard, strictAMIs]);
    ami_mean = mean([ami_standard, strictAMIs]);

    fprintf('\n============================================================\n');
    fprintf(' Strict validation summary\n');
    fprintf('============================================================\n');
    fprintf('    PAM AMI validated           = %.10f\n', ami_pam_val);
    fprintf('    Standard validated AMI      = %.10f\n', ami_standard);
    for k = 1:numel(strictResults)
        fprintf('    %-25s = %.10f\n', strictResults(k).name + " AMI", strictResults(k).ami);
    end
    fprintf('    Validation AMI range        = [%.10f, %.10f]\n', ami_min, ami_max);
    fprintf('    Validation AMI mean         = %.10f\n', ami_mean);
    fprintf('    Conservative gain min       = %.10f\n', ami_min - ami_pam_val);
    fprintf('    Optimistic gain max         = %.10f\n', ami_max - ami_pam_val);

    result = struct();
    result.params = params;
    result.x = x;
    result.p = p;
    result.ami_pam_val = ami_pam_val;
    result.ami_standard = ami_standard;
    result.fastTable = fastTable;
    result.strictResults = strictResults;
    result.ami_min = ami_min;
    result.ami_max = ami_max;
    result.ami_mean = ami_mean;

    save('result_validate_joint_solution_M8_SNR15_strict.mat', 'result');
    fprintf('\n[Save] Saved result_validate_joint_solution_M8_SNR15_strict.mat\n');

    fig = figure('Color','w','Name','Strict validation comparison');
    names = [{'standard'}, {strictResults.name}];
    vals = [ami_standard, strictAMIs];
    bar(vals);
    grid on; box on;
    set(gca, 'XTick', 1:numel(names), 'XTickLabel', names);
    ylabel('AMI [bits/symbol]');
    title('Strict validation comparison for final joint GS/PS solution');
    exportgraphics(fig, 'Strict_validation_joint_solution_M8_SNR15.png', 'Resolution', 300);
    fprintf('[Plot] Saved Strict_validation_joint_solution_M8_SNR15.png\n');

    fprintf('\nDone.\n');
end

% =========================================================================
% Local validation helpers
% =========================================================================
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

    % Add a denser positive mid region because IM/DD output is mostly positive.
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

function print_constraints(x, p, P_avg, minGap)
    fprintf('    sum(p)      = %.12f\n', sum(p));
    fprintf('    dot(p,x)    = %.12f\n', dot(p,x));
    fprintf('    target power= %.12f\n', P_avg);
    fprintf('    min(p)      = %.3e\n', min(p));
    fprintf('    min diff(x) = %.3e\n', min(diff(x)));
    fprintf('    minGap      = %.3e\n', minGap);
end
