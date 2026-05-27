%% run_GS_M8_SNR15.m
% Single-case geometric shaping simulation:
% M = 8, SNR = 15 dB
%
% Goal:
%   Maximize AMI using simulated annealing.
%   Report:
%     - Best GS constellation
%     - Validated AMI of Uniform PAM
%     - Validated AMI of GS
%     - AMI gain = AMI_GS - AMI_PAM

clear; clc; close all;
addpath(pwd);

fprintf('============================================================\n');
fprintf('   Geometric Shaping Simulation: M=8, SNR=15 dB\n');
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
cfg.SA.minGap       = 0.05;
cfg.SA.projectIters = 3;

% SA search settings
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

fprintf('[6] Running geometric shaping optimization...\n\n');

[out, results] = sa_multistart(cfg, x_pam(:), "[M8-SNR15]");

x_gs_fast = out.bestX(:).';
ami_gs_fast_best = out.bestMI;

fprintf('\n[6] Fast-SA result\n');
fprintf('    best start       = %d\n', out.bestStart);
fprintf('    GS x fast        = ['); fprintf(' %.6f', x_gs_fast); fprintf(' ]\n');
fprintf('    GS AMI fast      = %.8f bits/symbol\n\n', ami_gs_fast_best);


%% 7. Validate best GS constellation

fprintf('[7] Validating best GS constellation...\n');

% Project once more just to be safe
x_gs = AMI_functions.project_constellation_1D(x_gs_fast(:), cfg);
x_gs = x_gs(:).';

ami_gs_fast = AMI_functions.AMI_noCSI_fast_grid( ...
    x_gs, px, params, ghN_h, y_grid);

ami_gs_val = AMI_functions.AMI_noCSI_validate( ...
    x_gs, px, params);

gain_fast = ami_gs_fast - ami_pam_fast;
gain_val  = ami_gs_val  - ami_pam_val;

fprintf('    GS x final        = ['); fprintf(' %.6f', x_gs); fprintf(' ]\n');
fprintf('    mean(x_gs)        = %.12f\n', mean(x_gs));
fprintf('    min(x_gs)         = %.12f\n', min(x_gs));
fprintf('    min diff(x_gs)    = %.12f\n\n', min(diff(x_gs)));

fprintf('    GS AMI fast       = %.8f bits/symbol\n', ami_gs_fast);
fprintf('    GS AMI validated  = %.8f bits/symbol\n', ami_gs_val);
fprintf('    PAM AMI validated = %.8f bits/symbol\n', ami_pam_val);
fprintf('    AMI gain validated= %.8f bits/symbol\n\n', gain_val);


%% 8. Save results

result = struct();
result.params       = params;
result.cfg          = cfg;
result.x_pam        = x_pam;
result.x_gs         = x_gs;
result.ami_pam_fast = ami_pam_fast;
result.ami_pam_val  = ami_pam_val;
result.ami_gs_fast  = ami_gs_fast;
result.ami_gs_val   = ami_gs_val;
result.gain_fast    = gain_fast;
result.gain_val     = gain_val;
result.out          = out;
result.results      = results;

save('result_GS_M8_SNR15.mat', 'result');

fprintf('[8] Saved results to result_GS_M8_SNR15.mat\n\n');


%% 9. Optional plot

figure;
stem(1:params.M, x_pam, 'o', 'LineWidth', 1.5); hold on;
stem(1:params.M, x_gs, 'x', 'LineWidth', 1.5);
grid on; box on;
xlabel('Symbol index');
ylabel('Constellation level x');
title(sprintf('M=8, SNR=%g dB, \\sigma_R^2=%.2f', ...
    params.SNR_dB, params.sigma_X_sq));
legend('Uniform PAM', 'Geometric Shaping', 'Location', 'best');

fprintf('Done.\n');