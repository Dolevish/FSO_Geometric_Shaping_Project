function sim_shapingGain_vs_turbulence()
%   Fixed SNR, sweeps turbulence levels (sigma_X^2), optimizes 1D PAM via SA.
%
%   SA objective: AMI_functions.AMI_noCSI_fast_grid  (1D GH mixture + trapz on pre-built grid)
%     - Avoids triple-GH correlation bug that inflates MI at high SNR
%     - y-grid built ONCE per turbulence level, reused across all SA iters
%
%   Validation: AMI_functions.AMI_noCSI_validate  (adaptive integral via calculate_Py_given_x)
%     - Fully independent method (MATLAB's integral() — adaptive, high precision)

    clear; clc; close all;

    %% ====================================================================
    %  1. SIMULATION PARAMETERS
    % =====================================================================
    M       = 16;
    P_avg   = 1;
    SNR_dB  = 20;       

    % Turbulence levels (scintillation index)
    sigX_levels = [0, 0.1, 0.2, 0.3];
    level_names = ["No Turb (0)", "Very Weak (0.1)", "Weak (0.2)", "Weak-Mod (0.3)"];
    nLevels = numel(sigX_levels);

    % 1D GH order for h inside fast grid evaluator
    ghN_h = 40;

    % SA tuning
    sa_maxIter   = 50000;
    sa_nStarts   = 6;
    sa_useParfor = true;

    sigma_n_sq = P_avg / 10^(SNR_dB/10);

    fprintf('╔══════════════════════════════════════════════════════════════╗\n');
    fprintf('║   FSO Geometric Shaping  — Turbulence Sweep                 ║\n');
    fprintf('╠══════════════════════════════════════════════════════════════╣\n');
    fprintf('║  M=%d | P_avg=%.1f | SNR=%ddB | σ_n²=%.4f                 ║\n', M, P_avg, SNR_dB, sigma_n_sq);
    fprintf('║  Turbulence levels: %d | SA: %d starts × %d iter           ║\n', nLevels, sa_nStarts, sa_maxIter);
    fprintf('║  Evaluator: 1D-GH(N=%d) + trapz grid (no triple-GH)       ║\n', ghN_h);
    fprintf('╚══════════════════════════════════════════════════════════════╝\n\n');

    % --- Benchmark ---
    fprintf('Benchmarking fast evaluator...\n');
    cfg_bench = build_config(M, P_avg, SNR_dB, 0.5, ghN_h, sa_maxIter, sa_nStarts, false);
    x_bench = define_constellation(M, P_avg, true, "mean");
    Nrep = 50;
    t0b = tic;
    for rep = 1:Nrep
        cfg_bench.AMI_Evaluator(x_bench);
    end
    ms_per_call = toc(t0b) / Nrep * 1000;
    est_start = ms_per_call * sa_maxIter / 1000;
    fprintf('  %.2f ms/call → ~%.0fs per SA start\n', ms_per_call, est_start);
    fprintf('  Grid size: %d points\n\n', numel(cfg_bench.y_grid));

    %% ====================================================================
    %  2. PREALLOCATE
    % =====================================================================
    R = struct( ...
        'sigX',        num2cell(sigX_levels), ...
        'x_pam',       cell(1, nLevels), ...
        'x_opt',       cell(1, nLevels), ...
        'mi_pam_fast', num2cell(zeros(1,nLevels)), ...
        'mi_opt_fast', num2cell(zeros(1,nLevels)), ...
        'mi_pam_val',  num2cell(zeros(1,nLevels)), ...
        'mi_opt_val',  num2cell(zeros(1,nLevels)), ...
        'gain_bits',   num2cell(zeros(1,nLevels)), ...
        'rt_sa',       num2cell(zeros(1,nLevels)), ...
        'rt_val',      num2cell(zeros(1,nLevels)));

    t_total = tic;
    t_level_times = zeros(1, nLevels);

    %% ====================================================================
    %  3. MAIN SWEEP
    % =====================================================================
    for idx = 1:nLevels
        cur_sig = sigX_levels(idx);
        t_level = tic;

        eta_str = sweep_eta(idx, nLevels, t_level_times);
        fprintf('\n╔══════════════════════════════════════════════════════════╗\n');
        fprintf('║  [%d/%d]  σ_R² = %.2f  %-22s  %s\n', ...
            idx, nLevels, cur_sig, level_names(idx), eta_str);
        fprintf('╚══════════════════════════════════════════════════════════╝\n');
        fprintf('  Elapsed: %s\n', fmt_time(toc(t_total)));

        % --- Config (builds y-grid for this turbulence level) ---
        cfg = build_config(M, P_avg, SNR_dB, cur_sig, ...
            ghN_h, sa_maxIter, sa_nStarts, sa_useParfor);
        fprintf('  y-grid: %d pts [%.2f, %.2f]\n', ...
            numel(cfg.y_grid), cfg.y_grid(1), cfg.y_grid(end));

        % --- PAM baseline ---
        fprintf('\n  [1/3] PAM baseline...\n');
        x_pam = define_constellation(M, P_avg, true, "mean");
        mi_pam_fast = cfg.AMI_Evaluator(x_pam);
        fprintf('  ✦ PAM MI (grid-GH) = %.6f bits/sym\n', mi_pam_fast);

        % --- SA Optimization ---
        fprintf('\n  [2/3] SA (%d starts × %d iter)...\n', sa_nStarts, sa_maxIter);
        t_sa = tic;
        [out, sa_results] = sa_multistart(cfg, x_pam, sprintf("[sig%.1f]", cur_sig));
        rt_sa = toc(t_sa);

        x_opt       = out.bestX(:);
        mi_opt_fast = cfg.AMI_Evaluator(x_opt);

        % Per-start summary
        fprintf('\n  SA results:\n');
        fprintf('  %-6s %-12s %-8s\n', 'Start', 'MI(best)', 'Time');
        fprintf('  %s\n', repmat('-', 1, 28));
        for k = 1:numel(sa_results)
            marker = '';
            if sa_results(k).bestMI == out.bestMI, marker = ' ◄'; end
            fprintf('  %-6d %-12.6f %-8s%s\n', k, ...
                sa_results(k).bestMI, fmt_time(sa_results(k).runtime), marker);
        end
        fprintf('  ✦ OPT MI (grid-GH) = %.6f  |  gain = %+.4f\n', ...
            mi_opt_fast, mi_opt_fast - mi_pam_fast);
        fprintf('  ⏱  SA: %s\n', fmt_time(rt_sa));

        % Constellation
        fprintf('\n  Constellations:\n');
        fprintf('    PAM: '); fprintf('%7.4f ', sort(x_pam)); fprintf('\n');
        fprintf('    OPT: '); fprintf('%7.4f ', sort(x_opt)); fprintf('\n');
        gaps = diff(sort(x_opt));
        fprintf('    Gaps: '); fprintf('%6.4f ', gaps); fprintf('\n');
        if min(gaps) < 1e-2
            fprintf('  ⚠  min gap %.2e < 0.01 — alphabet reduction suspected\n', min(gaps));
        end

        % --- Validation (independent: adaptive integral + grid quadrature) ---
        fprintf('\n  [3/3] Validation (adaptive integral)...\n');
        t_val = tic;

        fprintf('    Grid quad PAM...');
        mi_pam_val = AMI_functions.AMI_noCSI_validate(x_pam, cfg.px, cfg);
        fprintf(' %.6f\n', mi_pam_val);

        fprintf('    Grid quad OPT...');
        mi_opt_val = AMI_functions.AMI_noCSI_validate(x_opt, cfg.px, cfg);
        fprintf(' %.6f\n', mi_opt_val);

        rt_val = toc(t_val);
        gain = mi_opt_val - mi_pam_val;

        % Cross-check fast vs validated
        d_pam = abs(mi_pam_fast - mi_pam_val);
        d_opt = abs(mi_opt_fast - mi_opt_val);

        fprintf('\n  ┌─ Cross-check ──────────────────────────────────┐\n');
        fprintf('  │  PAM: fast=%.5f  validated=%.5f  Δ=%.2e │\n', ...
            mi_pam_fast, mi_pam_val, d_pam);
        fprintf('  │  OPT: fast=%.5f  validated=%.5f  Δ=%.2e │\n', ...
            mi_opt_fast, mi_opt_val, d_opt);
        fprintf('  │  ★ SHAPING GAIN (validated) = %+.4f bits     │\n', gain);
        fprintf('  └───────────────────────────────────────────────┘\n');

        if max(d_pam, d_opt) > 0.05
            fprintf('  ⚠  Fast-vs-validated gap > 0.05 — increase ghN_h or grid density\n');
        end
        fprintf('  ⏱  Validation: %s\n', fmt_time(rt_val));

        % --- Store ---
        R(idx).x_pam       = x_pam(:);
        R(idx).x_opt       = x_opt(:);
        R(idx).mi_pam_fast = mi_pam_fast;
        R(idx).mi_opt_fast = mi_opt_fast;
        R(idx).mi_pam_val  = mi_pam_val;
        R(idx).mi_opt_val  = mi_opt_val;
        R(idx).gain_bits   = gain;
        R(idx).rt_sa       = rt_sa;
        R(idx).rt_val      = rt_val;

        t_level_times(idx) = toc(t_level);

        % Running progress
        fprintf('\n  ── Progress (%d/%d) ──\n', idx, nLevels);
        for ii = 1:idx
            fprintf('    σ²_R=%.1f: PAM=%.4f  OPT=%.4f  Gain=%+.4f\n', ...
                sigX_levels(ii), R(ii).mi_pam_val, R(ii).mi_opt_val, R(ii).gain_bits);
        end
        fprintf('  Level: %s  |  Total: %s', ...
            fmt_time(t_level_times(idx)), fmt_time(toc(t_total)));
        if idx < nLevels
            fprintf('  |  %s', sweep_eta(idx+1, nLevels, t_level_times));
        end
        fprintf('\n');
    end

    total_time = toc(t_total);

    %% ====================================================================
    %  4. SUMMARY TABLE
    % =====================================================================
    fprintf('\n\n');
    fprintf('╔═══════════════════════════════════════════════════════════════════════════════╗\n');
    fprintf('║                 RESULTS SUMMARY  (M=%d, SNR=%d dB)                  ║\n', M, SNR_dB);
    fprintf('╠═════════╦═══════════╦═══════════╦═══════════╦═══════════╦═════════╦══════════╣\n');
    fprintf('║  σ_R²   ║ PAM(fast) ║ OPT(fast) ║ PAM(val)  ║ OPT(val)  ║  Gain   ║ SA time  ║\n');
    fprintf('╠═════════╬═══════════╬═══════════╬═══════════╬═══════════╬═════════╬══════════╣\n');
    for idx = 1:nLevels
        fprintf('║  %5.2f  ║  %.4f   ║  %.4f   ║  %.4f   ║  %.4f   ║ %+.4f ║ %7s  ║\n', ...
            R(idx).sigX, R(idx).mi_pam_fast, R(idx).mi_opt_fast, ...
            R(idx).mi_pam_val, R(idx).mi_opt_val, R(idx).gain_bits, fmt_time(R(idx).rt_sa));
    end
    fprintf('╚═════════╩═══════════╩═══════════╩═══════════╩═══════════╩═════════╩══════════╝\n');
    fprintf('\nTotal: %s (%.1f min)\n', fmt_time(total_time), total_time/60);

    %% ====================================================================
    %  5. FIGURE 1: CONSTELLATION SUBPLOTS
    % =====================================================================
    nCols = min(3, nLevels);
    nRows = ceil(nLevels / nCols);

    fig1 = figure('Name', 'Constellation Comparison', 'Color', 'w', ...
        'Units', 'normalized', 'Position', [0.02 0.15 0.96 0.75]);

    for idx = 1:nLevels
        subplot(nRows, nCols, idx);
        hold on; grid on; box on;

        x_pam_s = sort(R(idx).x_pam);
        x_opt_s = sort(R(idx).x_opt);

        stem(x_pam_s, ones(size(x_pam_s)), 'bo', ...
            'MarkerSize', 8, 'LineWidth', 1.5, 'MarkerFaceColor', 'none');
        stem(x_opt_s, 0.5*ones(size(x_opt_s)), 'rx', ...
            'MarkerSize', 10, 'LineWidth', 2.0);

        ylim([-0.1, 1.4]);
        yticks([0.5, 1.0]); yticklabels({'Optimized', 'PAM'});
        xlabel('Intensity level');
        title(sprintf('\\sigma_R^2 = %.1f\nGain = %+.3f bits', ...
            R(idx).sigX, R(idx).gain_bits), 'FontSize', 10);

        xl = xlim;
        text(xl(2)*0.95, 1.3, sprintf('PAM: %.3f', R(idx).mi_pam_val), ...
            'FontSize', 8, 'HorizontalAlignment', 'right', 'Color', 'b');
        text(xl(2)*0.95, 1.15, sprintf('OPT: %.3f', R(idx).mi_opt_val), ...
            'FontSize', 8, 'HorizontalAlignment', 'right', 'Color', 'r');

        if idx == 1
            legend('Uniform PAM', 'GS Optimized', 'Location', 'northwest', 'FontSize', 7);
        end
        hold off;
    end
    sgtitle(sprintf('Geometric Shaping — M=%d, SNR=%d dB', M, SNR_dB), ...
        'FontSize', 13, 'FontWeight', 'bold');

    %% ====================================================================
    %  6. FIGURE 2: AMI + GAIN SUMMARY
    % =====================================================================
    fig2 = figure('Name', 'AMI Summary', 'Color', 'w', ...
        'Units', 'normalized', 'Position', [0.1 0.1 0.75 0.7]);

    subplot(1,2,1); hold on; grid on; box on;
    sigX_vec = [R.sigX]; pam_v = [R.mi_pam_val]; opt_v = [R.mi_opt_val];

    plot(sigX_vec, pam_v, 'b--o', 'LineWidth', 1.8, 'MarkerSize', 8, ...
        'MarkerFaceColor', 'b', 'DisplayName', 'Uniform PAM');
    plot(sigX_vec, opt_v, 'r-s', 'LineWidth', 2.0, 'MarkerSize', 9, ...
        'MarkerFaceColor', 'r', 'DisplayName', 'GS Optimized');
    yline(log2(M), 'k:', sprintf('log_2(%d)', M), ...
        'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

    xlabel('\sigma_R^2 (Rytov Variance)'); ylabel('AMI [bits/symbol]');
    title('AMI vs Turbulence'); legend('Location', 'southwest');
    xlim([min(sigX_vec)-0.05, max(sigX_vec)+0.1]);

    subplot(1,2,2); hold on; grid on; box on;
    gv = [R.gain_bits];
    b = bar(1:nLevels, gv, 0.6); b.FaceColor = 'flat';
    for k = 1:nLevels
        if gv(k)>0.01, b.CData(k,:)=[0.85 0.33 0.10];
        else,          b.CData(k,:)=[0.5 0.5 0.5]; end
    end
    xtl = arrayfun(@(v)sprintf('%.1f',v), sigX_levels, 'Uni', false);
    xticks(1:nLevels); xticklabels(xtl);
    xlabel('\sigma_R^2'); ylabel('Shaping Gain [bits/sym]');
    title('Shaping Gain (OPT - PAM)');
    for k=1:nLevels
        text(k, gv(k)+0.005, sprintf('%+.3f',gv(k)), ...
            'HorizontalAlignment','center','FontSize',9,'FontWeight','bold');
    end
    sgtitle(sprintf('GS Summary — M=%d, SNR=%d dB', M, SNR_dB), ...
        'FontSize', 13, 'FontWeight', 'bold');

    %% ====================================================================
    %  7. FIGURE 3: RUNTIME
    % =====================================================================
    fig3 = figure('Name', 'Runtime', 'Color', 'w', ...
        'Units', 'normalized', 'Position', [0.25 0.2 0.5 0.5]);
    bh = bar(1:nLevels, [[R.rt_sa];[R.rt_val]]', 'grouped');
    bh(1).FaceColor = [0 0.45 0.74]; bh(2).FaceColor = [0.47 0.67 0.19];
    xticks(1:nLevels); xticklabels(xtl);
    xlabel('\sigma_R^2'); ylabel('Time [s]');
    title(sprintf('Runtime (Total: %s)', fmt_time(total_time)));
    legend('SA','Validation','Location','northwest');
    grid on; box on;

    fprintf('\n✓ Done.\n');
end


% =========================================================================
%  BUILD CONFIG
% =========================================================================
function cfg = build_config(M, P_avg, SNR_dB, sigma_X_sq, ...
        ghN_h, sa_maxIter, sa_nStarts, sa_useParfor)

    cfg = struct();
    cfg.M          = M;
    cfg.P_avg      = P_avg;
    cfg.SNR_dB     = SNR_dB;
    cfg.sigma_X_sq = sigma_X_sq;
    cfg.sigma_n_sq = P_avg / (10^(SNR_dB / 10));
    cfg.R          = 1;

    % Log-normal fading (scintillation index convention)
    cfg.sig_t = sqrt(log(1 + sigma_X_sq));
    cfg.mu_t  = -0.5 * cfg.sig_t^2;

    cfg.px = ones(M, 1) / M;

    % --- Pre-build y-grid for this turbulence level ---
    % Use generous x_max bound 
    x_max_bound = 5;   % generous for mean(x)=1, M=8
    cfg.ghN_h  = ghN_h;
    cfg.y_grid = AMI_functions.build_noCSI_y_grid(cfg, x_max_bound);

    % SA
    sa = struct();
    sa.maxIter      = sa_maxIter;
    sa.itersPerTemp = 50;
    sa.T0           = 0.4;
    sa.Tf           = 1e-3;
    sa.nBlocks      = ceil(sa.maxIter / sa.itersPerTemp);
    sa.coolingRate   = exp(log(sa.Tf/sa.T0) / sa.nBlocks);

    sa.baseStd0 = 0.15; sa.baseStdMin = 1e-3; sa.baseStdMax = 0.5;
    sa.targetAccLo = 0.20; sa.targetAccHi = 0.60;
    sa.baseStdGrow = 1.25; sa.baseStdShrink = 0.80;

    sa.enforce_sort = true; sa.imdd_mode = true;
    sa.powerConstraint = "mean"; sa.enforce_power = true;
    sa.minGap = 0.05; sa.projectIters = 2;

    sa.nStarts = sa_nStarts; sa.useParallel = sa_useParfor;
    sa.numWorkers = []; sa.closePoolWhenDone = false;

    % Progress logging: every 5000 iter 
    sa.logEvery = 5000;

    cfg.SA = sa;

    % --- Objective: grid-based evaluator ---
    y_grid_local = cfg.y_grid;
    ghN_local    = ghN_h;
    cfg.AMI_Evaluator = @(x_in) eval_fast(x_in, cfg, ghN_local, y_grid_local);
end


function mi = eval_fast(x_in, cfg, ghN_h, y_grid)
    x  = AMI_functions.project_constellation_1D(x_in, cfg);
    mi = AMI_functions.AMI_noCSI_fast_grid(x, cfg.px, cfg, ghN_h, y_grid);
end


% =========================================================================
%  HELPERS
% =========================================================================
function s = fmt_time(sec)
    if sec < 60,       s = sprintf('%.1fs', sec);
    elseif sec < 3600, s = sprintf('%dm%02.0fs', floor(sec/60), mod(sec,60));
    else,              s = sprintf('%dh%02dm', floor(sec/3600), floor(mod(sec,3600)/60));
    end
end

function s = sweep_eta(idx, nLevels, times)
    done = times(1:idx-1); done = done(done>0);
    if isempty(done), s = 'ETA --'; return; end
    s = sprintf('ETA ~%s', fmt_time(mean(done)*(nLevels-idx+1)));
end