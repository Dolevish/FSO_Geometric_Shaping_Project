% main_validation_channel.m
%
% PURPOSE: Channel-level validation for the FSO no-CSI geometric shaping project.
%   Verifies P(y|x), APPs, decision boundaries, SER and 
%   the two AMI evaluators actually used in production:
%       • AMI_functions.AMI_noCSI_fast_grid  ← SA objective (what gets optimised)
%       • AMI_functions.AMI_noCSI_validate   ← ground-truth validator
%   A large Δ between the two flags evaluator inaccuracy
%
%   Dependencies:
%       - AMI_functions.m
%       - BER_functions.m
%       - calculate_Py_given_x.m


%% ── 1. PARAMETERS ──────────────────────────────────────────────────────────
clear; clc; close all;
addpath(pwd);

params.R          = 1;      % Detector responsivity [A/W]
params.sigma_X_sq = 0.1;    % Rytov variance (weak turbulence parameter)
params.M          = 4;      % Constellation size
params.P_avg      = 1;      % Average optical power [normalised]
params.SNR_dB     = 15;     % Target SNR [dB]
params.sigma_n_sq = params.P_avg / 10^(params.SNR_dB/10);

% Log-normal fading parameters  
sig_t_sq       = log(1 + params.sigma_X_sq);
params.sig_t   = sqrt(sig_t_sq);
params.mu_t    = -0.5 * sig_t_sq;

fprintf('╔══════════════════════════════════════════════════════════════╗\n');
fprintf('║         FSO Channel Validation Script                        ║\n');
fprintf('╠══════════════════════════════════════════════════════════════╣\n');
fprintf('║  M=%-2d | P_avg=%.1f | SNR=%ddB | σ_n²=%.4f             ║\n', ...
    params.M, params.P_avg, params.SNR_dB, params.sigma_n_sq);
fprintf('║  σ_R²=%.2f | μ_t=%.4f | σ_t=%.4f                     ║\n', ...
    params.sigma_X_sq, params.mu_t, params.sig_t);
fprintf('╚══════════════════════════════════════════════════════════════╝\n\n');

if ~exist('calculate_Py_given_x','file')
    error('calculate_Py_given_x.m not found in path.');
end


%% ── 2. IM/DD CONSTELLATION ──────────────────────────────────────────────────
%   Non-negative uniform PAM — identical to what SA starts from.
%   x(1)=0 pinned at origin (IM/DD boundary condition).
fprintf('[2] Building IM/DD PAM constellation...\n');

step = 2*params.P_avg / (params.M - 1);   % spacing for mean = P_avg
constellation_points = (0:params.M-1) * step;

assert(abs(mean(constellation_points) - params.P_avg) < 1e-10, ...
    'Mean power constraint violated.');

symbol_probs = ones(1, params.M) / params.M;

fprintf('  Constellation (IM/DD PAM-%d):\n', params.M);
fprintf('    x = [');  fprintf(' %.4f', constellation_points);  fprintf(' ]\n');
fprintf('  Mean power: %.4f  (target: %.4f)\n\n', ...
    mean(constellation_points), params.P_avg);


%% ── 3. Y-GRIDS ──────────────────────────────────────────────────────────────
%   (A) PLOT grid  — uniform, moderate width, for clean figures
%   (B) SA grid    — non-uniform, from AMI_functions.build_noCSI_y_grid

fprintf('[3] Building y-grids...\n');

noise_std    = sqrt(params.sigma_n_sq);
h_hi_plot    = logninv(1-1e-3, params.mu_t, params.sig_t);
y_abs_plot   = params.R * max(constellation_points) * h_hi_plot + 7*noise_std;
N_plot       = 2001;
y_range_plot = linspace(-0.5, y_abs_plot * 1.05, N_plot);

x_max_bound = max(constellation_points) * 1.5;
y_sa_grid   = AMI_functions.build_noCSI_y_grid(params, x_max_bound);

fprintf('  (A) Plot grid  (uniform):     [%.2f, %.2f]  N=%d\n', ...
    y_range_plot(1), y_range_plot(end), N_plot);
fprintf('  (B) SA grid    (non-uniform): [%.2f, %.2f]  N=%d\n\n', ...
    y_sa_grid(1), y_sa_grid(end), numel(y_sa_grid));


%% ── 4. P(y|x) ON PLOT GRID ─────────────────────────────────────────────────
fprintf('[4] Computing P(y|x) on plot grid... ');

Py_given_x_plot = zeros(params.M, N_plot);
for i = 1:params.M
    Py_given_x_plot(i,:) = calculate_Py_given_x( ...
        y_range_plot, constellation_points(i), params);
end
Py_total_plot = symbol_probs * Py_given_x_plot;

fprintf('Done.\n\n');


%% ── 5. FIGURE 1: P(y|x) AND P(y) ───────────────────────────────────────────
fprintf('[5] Plotting P(y|x) and P(y)...\n');

figure(1); clf;
hold on;
clrs = lines(params.M);
leg  = cell(1, params.M+1);

for i = 1:params.M
    plot(y_range_plot, Py_given_x_plot(i,:), ...
        '--', 'Color', clrs(i,:), 'LineWidth', 1.5);
    leg{i} = sprintf('P(y|x=%.2f)', constellation_points(i));
end

% P(y) total 
plot(y_range_plot, Py_total_plot, '-', 'Color', [0 0 0], 'LineWidth', 2.5);
leg{params.M+1} = 'P(y)  [Total]';

hold off; grid on; box on;
xlabel('Received signal y');
ylabel('Probability density');
title(sprintf('P(y|x) — IM/DD PAM-%d  (SNR=%ddB, \\sigma_R^2=%.2f)', ...
    params.M, params.SNR_dB, params.sigma_X_sq), 'FontSize', 12);
legend(leg, 'Location', 'best', 'FontSize', 9);
ylim([0, max(Py_given_x_plot(:)) * 1.2]);


%% ── 6. APPs: P(xi|y) ────────────────────────────────────────────────────────
fprintf('[6] Computing APPs P(xi|y)  (log-domain)... ');

logPyx   = log(max(Py_given_x_plot, realmin));
logPrior = log(symbol_probs(:));
logPost  = bsxfun(@plus, logPyx, logPrior);
maxL     = max(logPost, [], 1);
logDen   = maxL + log(sum(exp(bsxfun(@minus, logPost, maxL)), 1));
P_xi_given_y_plot = exp(bsxfun(@minus, logPost, logDen));   % M × N_plot

fprintf('Done.\n\n');


%% ── 7. DECISION BOUNDARIES (using BER_functions) ───────────────────────────
fprintf('[7] Finding and refining decision boundaries...\n');

boundaries = -inf;
for i = 1:(params.M-1)
    y_mid = params.R * (constellation_points(i) + constellation_points(i+1)) / 2;
    if params.sigma_X_sq <= 1e-12
        y_b = y_mid;
    else
        try
            % Use BER_functions for grid-based search
            y_b = BER_functions.find_boundary_loglike(y_range_plot, ...
                Py_given_x_plot(i,:), Py_given_x_plot(i+1,:), y_mid);
        catch ME
            warning('Boundary search failed (%s). Using midpoint.', ME.message);
            y_b = y_mid;
        end
    end
    boundaries(end+1) = y_b; 
end
boundaries(end+1) = inf;

tolB       = max(1e-9, 0.001 * min(diff(constellation_points)));
boundaries = sort(unique(boundaries));
boundaries = boundaries([true, diff(boundaries) > tolB]);

% Refine using BER_functions
for i = 1:(params.M-1)
    [y_ref, ~] = BER_functions.refine_boundary_continuous(boundaries(i+1), ...
        constellation_points(i), constellation_points(i+1), ...
        constellation_points, symbol_probs, params);
    boundaries(i+1) = y_ref;
end

fprintf('  Refined decision boundaries:\n');
for i = 2:numel(boundaries)-1
    fprintf('    between x=%.4f and x=%.4f  →  y* = %.6f\n', ...
        constellation_points(i-1), constellation_points(i), boundaries(i));
end
fprintf('\n');


%% ── 8. FIGURE 2: APPs WITH DECISION REGIONS ─────────────────────────────────
fprintf('[8] Plotting APPs with decision region boundaries...\n');

figure(2); clf;
hold on;
leg_app = cell(1, params.M);
for i = 1:params.M
    plot(y_range_plot, P_xi_given_y_plot(i,:), ...
        '-', 'Color', clrs(i,:), 'LineWidth', 1.8);
    leg_app{i} = sprintf('P(x=%.2f | y)', constellation_points(i));
end

yl = [0, 1.05];
for k = 2:(numel(boundaries)-1)
    if isfinite(boundaries(k))
        plot([boundaries(k) boundaries(k)], yl, 'k--', 'LineWidth', 1.2, ...
            'HandleVisibility', 'off');
    end
end

% Annotation explaining the dashed lines
ann_x = boundaries(2) + 0.03*(y_range_plot(end)-y_range_plot(1));
text(ann_x, 0.87, {'Dashed black lines', '= decision region', 'boundaries'}, ...
    'FontSize', 8.5, 'Color', [0.15 0.15 0.15], ...
    'BackgroundColor', 'w', 'EdgeColor', [0.6 0.6 0.6]);

hold off; grid on; box on;
ylim(yl);
xlabel('Received signal y');
ylabel('A posteriori probability P(xi|y)');
title(sprintf('APPs — IM/DD PAM-%d  (SNR=%ddB, \\sigma_R^2=%.2f)', ...
    params.M, params.SNR_dB, params.sigma_X_sq), 'FontSize', 12);
legend(leg_app, 'Location', 'east', 'FontSize', 9);


%% ── 9. SYMBOL ERROR RATE (using BER_functions) ──────────────────────────────
fprintf('[9] Computing SER (quadgk over refined decision regions)...\n');

[avg_ser, symbol_error_probs] = BER_functions.compute_SER_quadgk( ...
    constellation_points, symbol_probs, boundaries, params);

fprintf('  Per-symbol P(e|x):\n');
for i = 1:params.M
    fprintf('    x=%.4f : %.4e\n', constellation_points(i), symbol_error_probs(i));
end
fprintf('  Average SER: %.4e\n\n', avg_ser);

% Also compute BER using the main function for comparison
ber_ml = BER_functions.calculate_BER_noCSI_ML(constellation_points(:), params);
fprintf('  BER (ML thresholds, GH integration): %.4e\n\n', ber_ml);


%% ── 10. AMI — PRODUCTION EVALUATORS ─────────────────────────────────────────
%
%  Uses EXACTLY the two evaluators that matter:
%    (a) AMI_noCSI_fast_grid  — SA objective (1D-GH mixture + trapz on SA grid)
%    (b) AMI_noCSI_validate   — ground truth (adaptive integral via calculate_Py_given_x)
%
%  How to interpret Δ(fast − validated):
%    < 0.005  →  evaluators agree; 
%    0.005–0.05  →  moderate bias; 
%    > 0.05   →  large bias; 
%
fprintf('[10] Computing AMI (production evaluators)...\n');

px    = symbol_probs(:);
ghN_h = 40;   % same as SA default

mi_fast = AMI_functions.AMI_noCSI_fast_grid( ...
    constellation_points, px, params, ghN_h, y_sa_grid);

mi_val = AMI_functions.AMI_noCSI_validate( ...
    constellation_points, px, params);

delta = abs(mi_fast - mi_val);
Hmax  = log2(params.M);

if     delta < 0.005, flag = 'Evaluators agree';
elseif delta < 0.05,  flag = 'Moderate gap — monitor at high SNR';
else,                 flag = 'Large gap';
end

fprintf('\n');
fprintf('  ┌─ AMI Results ──────────────────────────────────────────┐\n');
fprintf('  │  (a) Fast  [SA objective]   : %.6f bits/sym          │\n', mi_fast);
fprintf('  │  (b) Validated [ground tth] : %.6f bits/sym          │\n', mi_val);
fprintf('  │                                                         │\n');
fprintf('  │  Δ (fast − validated)       : %.2e bits              │\n', delta);
fprintf('  │  H_max = log2(%d)            : %.6f bits/sym          │\n', params.M, Hmax);
fprintf('  │  %s\n', [flag, repmat(' ', 1, max(0, 50-length(flag))), '│']);
fprintf('  └─────────────────────────────────────────────────────────┘\n\n');


%% ── 11. VALIDATION TESTS ────────────────────────────────────────────────────
fprintf('[11] Validation tests\n');
fprintf('══════════════════════════════════════════════════════════════\n');

% T1: APP columns sum to 1
err_sum = max(abs(sum(P_xi_given_y_plot, 1) - 1));
fprintf('[T1] APP normalisation  max|Σ_x P(x|y) - 1|\n');
fprintf('     max error = %.3e    %s  (threshold: 1e-10)\n\n', ...
    err_sum, pass_fail(err_sum < 1e-10));

% T2: APP equality at each decision boundary
%   Evaluate P(y|x) DIRECTLY at exact y* via continuous channel law —
%   not via grid lookup (dsearchn would introduce discretisation error).
fprintf('[T2] APP equality at decision boundaries\n');
all_ok = true;
for i = 1:(params.M-1)
    yb = boundaries(i+1);
    if ~isfinite(yb), continue; end
    p_all = arrayfun(@(xx) max(calculate_Py_given_x(yb, xx, params), realmin), ...
                     constellation_points(:)');
    post  = (symbol_probs .* p_all) / sum(symbol_probs .* p_all);
    app_i = post(i);  app_j = post(i+1);
    gap   = abs(app_i - app_j);
    ok    = gap < 1e-3;
    all_ok = all_ok && ok;
    fprintf('     y*=%.4f : P(x=%.3f|y*)=%.4f  P(x=%.3f|y*)=%.4f  |Δ|=%.2e  %s\n', ...
        yb, constellation_points(i), app_i, constellation_points(i+1), app_j, gap, pass_fail(ok));
end
fprintf('     Overall: %s\n\n', pass_fail(all_ok));

% T3: Boundaries strictly increasing
fin_b = boundaries(~isinf(boundaries));
diffs = diff(fin_b);
fprintf('[T3] Decision boundaries strictly increasing\n');
fprintf('     [');  fprintf(' %.4f', fin_b);  fprintf(' ]\n');
fprintf('     Min gap = %.4e    %s\n\n', min(diffs), pass_fail(all(diffs > 0)));

% T5: AMI in valid range
in_range = (mi_val >= -1e-9) && (mi_val <= Hmax + 1e-9);
fprintf('[T5] AMI(validated) ∈ [0, log2(M)]\n');
fprintf('     %.6f ∈ [0, %.4f]    %s\n\n', mi_val, Hmax, pass_fail(in_range));

% Hard asserts
assert(err_sum  < 1e-10, 'APP columns must sum to 1');
assert(in_range,         'AMI out of [0, log2(M)]');
assert(all(diffs > 0),   'Decision boundaries must be strictly increasing');

fprintf('══════════════════════════════════════════════════════════════\n');
fprintf('✓ All hard assertions passed.\n');


%% ═══════════════════════════════════════════════════════════════════════════
%  HELPER FUNCTIONS
% ════════════════════════════════════════════════════════════════════════════

function s = pass_fail(cond)
    if cond, s = '✓ PASS'; else, s = '✗ FAIL'; end
end
