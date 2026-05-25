function sim_AMI_vs_SNR()
%SIM_AMI_VS_SNR  AMI vs SNR for FSO geometric shaping.
%
%   Grid sweep: SNR × turbulence levels (σ²_R).
%   For each point: SA optimizes 1D PAM constellation, validated AMI computed.
%   Parallel over grid points (parfor) — NO nested parfor; SA restarts serial.
%
%   Evaluator:  AMI_functions.AMI_noCSI_fast_grid  (1D-GH + trapz grid)
%   Validation: AMI_functions.AMI_noCSI_validate   (adaptive integral)
%
%   Output:
%     - Figure 1: AMI vs SNR curves — PAM (dashed) vs GS-optimized (solid)
%     - Figure 2: Shaping Gain vs SNR (bits/symbol)
%     - Figure 3: BER vs SNR (semilog scale) — No-CSI receiver, ML thresholds
%     - Figure 4: Shaping Gain vs SNR (dB)
%
%   Dependencies:
%     - AMI_functions.m
%     - BER_functions.m
%     - calculate_Py_given_x.m
%     - define_constellation.m
%     - simulated_annealing.m

    clear; clc; close all;

    %% ====================================================================
    %  1. SIMULATION PARAMETERS
    % =====================================================================
    M     = 4;
    P_avg = 1;

    snr_vec        = 5:5:30;
    sigX_vec       = [0.0,  0.1,  0.2,  0.3];
    level_names    = ["No Turb (0)", "Very Weak (0.1)", "Weak (0.2)", "Weak-Mod (0.3)"];

    nSNR    = numel(snr_vec);
    nLevels = numel(sigX_vec);
    nPoints = nSNR * nLevels;

    ghN_h      = 40;    % GH nodes for h integral (1D) 
    sa_maxIter = 10000; % per start
    sa_nStarts = 8;     % serial restarts per grid point 

    fprintf('╔══════════════════════════════════════════════════════════════╗\n');
    fprintf('║   FSO Geometric Shaping — AMI vs SNR Sweep                   ║\n');
    fprintf('╠══════════════════════════════════════════════════════════════╣\n');
    fprintf('║  M=%d | P_avg=%.1f | SNR: %d→%d dB (%d pts)                 ║\n', ...
        M, P_avg, snr_vec(1), snr_vec(end), nSNR);
    fprintf('║  Turbulence levels: %d | Grid points: %d                      ║\n', nLevels, nPoints);
    fprintf('║  SA: %d starts × %d iter (per point) | GH N=%d             ║\n', ...
        sa_nStarts, sa_maxIter, ghN_h);
    fprintf('║  Evaluator: 1D-GH + trapz grid (no triple-GH)              ║\n');
    fprintf('║  BER: No-CSI receiver with ML-optimal thresholds           ║\n');
    fprintf('╚══════════════════════════════════════════════════════════════╝\n\n');

    %% ====================================================================
    %  2. BENCHMARK
    % =====================================================================
    fprintf('Benchmarking fast evaluator (SNR=20dB, σ²=0.1)...\n');
    cfg_bench = build_config(M, P_avg, 20, 0.1, ghN_h, sa_maxIter, 1, false);
    x_bench   = define_constellation(M, P_avg, true, "mean");
    Nrep = 30; t0b = tic;
    for rep = 1:Nrep, cfg_bench.AMI_Evaluator(x_bench); end
    ms_per_call = toc(t0b)/Nrep*1000;
    est_pt = ms_per_call * sa_maxIter * sa_nStarts / 1000;
    fprintf('  %.2f ms/call → ~%s per grid point (%d starts × %d iter)\n', ...
        ms_per_call, fmt_time(est_pt), sa_nStarts, sa_maxIter);
    fprintf('  Grid size: %d points\n\n', numel(cfg_bench.y_grid));

    %% ====================================================================
    %  3. FLATTEN (SNR, σ²) GRID FOR PARFOR
    %     Row-major: task k = (i_sig-1)*nSNR + j_snr
    %     → reshape(result, nLevels, nSNR) gives rows=turbulence, cols=SNR
    % =====================================================================
    [SNR_Grid, SIG_Grid] = meshgrid(snr_vec, sigX_vec);
    task_snr = SNR_Grid(:);
    task_sig = SIG_Grid(:);

    %% ====================================================================
    %  4. PARALLEL POOL
    % =====================================================================
    pool = gcp('nocreate');
    if isempty(pool)
        pool = parpool('Processes');
    end
    nWorkers = pool.NumWorkers;
    fprintf('Parallel pool: %d workers\n', nWorkers);
    est_total = est_pt * ceil(nPoints / nWorkers);
    fprintf('Estimated total: ~%s\n\n', fmt_time(est_total));

    %% ====================================================================
    %  5. PREALLOCATE
    % =====================================================================
    res_pam_fast = zeros(nPoints, 1);
    res_opt_fast = zeros(nPoints, 1);
    res_pam_val  = zeros(nPoints, 1);
    res_opt_val  = zeros(nPoints, 1);
    res_rt_sa    = zeros(nPoints, 1);
    
    % BER results
    res_ber_pam  = zeros(nPoints, 1);
    res_ber_opt  = zeros(nPoints, 1);
    
    % Store optimal constellations for reproducibility
    res_x_opt    = cell(nPoints, 1);

    %% ====================================================================
    %  6. PARFOR OVER GRID POINTS
    % =====================================================================
    Q = parallel.pool.DataQueue;
    count_done = 0;
    afterEach(Q, @(~) progress_tick());
    function progress_tick()
        count_done = count_done + 1;
        fprintf('  ▸ %2d/%d done\n', count_done, nPoints);
    end

    fprintf('Running grid (%d×%d = %d points)...\n', nLevels, nSNR, nPoints);
    t_total = tic;

    parfor k = 1:nPoints
        cur_snr = task_snr(k);
        cur_sig = task_sig(k);

        % Config with silent SA (logEvery=0 avoids interleaved parfor output)
        cfg_k = build_config(M, P_avg, cur_snr, cur_sig, ghN_h, sa_maxIter, 1, false);
        cfg_k.SA.logEvery = 0;

        % PAM baseline
        x_pam_k = define_constellation(M, P_avg, true, "mean");
        mi_pam_fast_k = cfg_k.AMI_Evaluator(x_pam_k);

        % SA — serial restarts (no nested parfor allowed inside parfor)
        t_sa = tic;
        best_mi = -inf;
        best_x  = x_pam_k(:);
        label_k = sprintf('[SNR%02d|s%.1f]', cur_snr, cur_sig);

        % Start 1: PAM init
        out1 = simulated_annealing(cfg_k, x_pam_k, 1, label_k);
        if out1.bestMI > best_mi
            best_mi = out1.bestMI;  best_x = out1.x_best;
        end

        % Starts 2..nStarts: random perturbations
        for s = 2:sa_nStarts
            x0_s = x_pam_k(:) + 0.15 * randn(M, 1);
            out_s = simulated_annealing(cfg_k, x0_s, s, label_k);
            if out_s.bestMI > best_mi
                best_mi = out_s.bestMI;  best_x = out_s.x_best;
            end
        end
        rt_sa_k = toc(t_sa);

        mi_opt_fast_k = cfg_k.AMI_Evaluator(best_x);

        % Validation (independent method — adaptive integral)
        mi_pam_val_k = AMI_functions.AMI_noCSI_validate(x_pam_k, cfg_k.px, cfg_k);
        mi_opt_val_k = AMI_functions.AMI_noCSI_validate(best_x,  cfg_k.px, cfg_k);

        % BER calculation (No-CSI with ML-optimal thresholds)
        ber_pam_k = BER_functions.calculate_BER_noCSI_ML(x_pam_k(:), cfg_k);
        ber_opt_k = BER_functions.calculate_BER_noCSI_ML(best_x(:),  cfg_k);

        % Store
        res_pam_fast(k) = mi_pam_fast_k;
        res_opt_fast(k) = mi_opt_fast_k;
        res_pam_val(k)  = mi_pam_val_k;
        res_opt_val(k)  = mi_opt_val_k;
        res_rt_sa(k)    = rt_sa_k;
        res_ber_pam(k)  = ber_pam_k;
        res_ber_opt(k)  = ber_opt_k;
        res_x_opt{k}    = best_x(:);

        send(Q, 1);
    end

    total_time = toc(t_total);
    fprintf('\n✓ Grid complete. Total: %s (%.1f min)\n\n', ...
        fmt_time(total_time), total_time/60);

    %% ====================================================================
    %  7. RESHAPE → (nLevels × nSNR)
    % =====================================================================
    ami_pam_fast = reshape(res_pam_fast, nLevels, nSNR);
    ami_opt_fast = reshape(res_opt_fast, nLevels, nSNR);
    ami_pam_val  = reshape(res_pam_val,  nLevels, nSNR);
    ami_opt_val  = reshape(res_opt_val,  nLevels, nSNR);
    ber_pam      = reshape(res_ber_pam,  nLevels, nSNR);
    ber_opt      = reshape(res_ber_opt,  nLevels, nSNR);
    
    % Reshape constellations to cell array (nLevels × nSNR)
    x_opt_grid   = reshape(res_x_opt, nLevels, nSNR);

    %% ====================================================================
    %  8. SUMMARY TABLES (one per turbulence level)
    % =====================================================================
    fprintf('\n\n');
    for i = 1:nLevels
        fprintf('╔═══════════════════════════════════════════════════════════════════════════════════════════════╗\n');
        fprintf('║  σ_R² = %.2f  —  %-38s                                    ║\n', sigX_vec(i), level_names(i));
        fprintf('╠═════════╦══════════╦══════════╦══════════╦══════════╦══════════╦══════════╦════════╦═════════╣\n');
        fprintf('║  SNR    ║ PAM(fast)║ OPT(fast)║ PAM(val) ║ OPT(val) ║  Gain    ║ Gain(dB) ║ Δ fast ║ BER imp ║\n');
        fprintf('╠═════════╬══════════╬══════════╬══════════╬══════════╬══════════╬══════════╬════════╬═════════╣\n');
        for j = 1:nSNR
            gain_bits = ami_opt_val(i,j) - ami_pam_val(i,j);
            gain_dB   = 10*log10(ami_opt_val(i,j) / max(ami_pam_val(i,j), 1e-6));
            d_fast    = abs(ami_opt_fast(i,j) - ami_opt_val(i,j));
            ber_ratio = ber_pam(i,j) / max(ber_opt(i,j), 1e-15);
            fprintf('║  %5.1fdB ║  %.4f  ║  %.4f  ║  %.4f  ║  %.4f  ║ %+.4f  ║ %+6.2fdB ║ %.4f ║  x%.1f   ║\n', ...
                snr_vec(j), ami_pam_fast(i,j), ami_opt_fast(i,j), ...
                ami_pam_val(i,j), ami_opt_val(i,j), gain_bits, gain_dB, d_fast, ber_ratio);
        end
        fprintf('╚═════════╩══════════╩══════════╩══════════╩══════════╩══════════╩══════════╩════════╩═════════╝\n\n');
    end

    % Shaping gain matrix summary
    fprintf('╔══════════════════════════════════════════════════════════════╗\n');
    fprintf('║            SHAPING GAIN (bits/sym)  [OPT(val) − PAM(val)]  ║\n');
    fprintf('╠══════════════╦');
    fprintf('%s', repmat('═════════╦', 1, nSNR-1)); fprintf('═════════╣\n');
    fprintf('║  σ²_R \\ SNR   ║');
    for j = 1:nSNR, fprintf(' %5.0fdB  ║', snr_vec(j)); end; fprintf('\n');
    fprintf('╠══════════════╬');
    fprintf('%s', repmat('═════════╬', 1, nSNR-1)); fprintf('═════════╣\n');
    for i = 1:nLevels
        fprintf('║  %.2f         ║', sigX_vec(i));
        for j = 1:nSNR
            fprintf('  %+.4f ║', ami_opt_val(i,j) - ami_pam_val(i,j));
        end
        fprintf('\n');
    end
    fprintf('╚══════════════╩');
    fprintf('%s', repmat('═════════╩', 1, nSNR-1)); fprintf('═════════╝\n');
    fprintf('\nTotal: %s (%.1f min)\n', fmt_time(total_time), total_time/60);

    %% ====================================================================
    %  9. COLORS AND MARKERS
    % =====================================================================
    colors = {
        [0.00  0.45  0.74],   % blue
        [0.64  0.08  0.18],   % dark red
        [0.47  0.67  0.19],   % green
        [0.93  0.69  0.13]    % gold/orange
    };
    markers = {'s', '*', 'o', 'd'};

    %% ====================================================================
    %  10. FIGURE 1: AMI vs SNR  
    % =====================================================================
    fig1 = figure('Name', 'AMI vs SNR', 'Color', 'w', ...
        'Units', 'normalized', 'Position', [0.02 0.52 0.45 0.40]);
    hold on; grid on; box on;

    for i = 1:nLevels
        c = colors{i};  m = markers{i};
        plot(snr_vec, ami_opt_val(i,:), ['-' m], ...
            'Color', c, 'LineWidth', 2.0, ...
            'MarkerFaceColor', c, 'MarkerSize', 7, ...
            'DisplayName', sprintf('GS opt  (σ²_R=%.2f)', sigX_vec(i)));
        plot(snr_vec, ami_pam_val(i,:), ['--' m], ...
            'Color', c, 'LineWidth', 1.5, ...
            'MarkerFaceColor', 'w', 'MarkerSize', 6, ...
            'DisplayName', sprintf('PAM    (σ²_R=%.2f)', sigX_vec(i)));
    end

    yline(log2(M), 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
    text(snr_vec(end)+0.3, log2(M), sprintf('log_2(%d)=%d', M, log2(M)), ...
        'FontSize', 9, 'VerticalAlignment', 'bottom');

    xlabel('SNR (dB)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('AMI (bits/symbol)', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('FSO Geometric Shaping — M=%d PAM', M), ...
        'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'southeast', 'NumColumns', 2, 'FontSize', 10);
    ylim([0, log2(M) * 1.08]);
    xlim([snr_vec(1)-1, snr_vec(end)+1]);
    xticks(snr_vec);
    hold off;

    %% ====================================================================
    %  11. FIGURE 2: SHAPING GAIN vs SNR (bits/symbol)
    % =====================================================================
    fig2 = figure('Name', 'Shaping Gain vs SNR', 'Color', 'w', ...
        'Units', 'normalized', 'Position', [0.50 0.52 0.45 0.40]);
    hold on; grid on; box on;

    for i = 1:nLevels
        if sigX_vec(i) == 0, continue; end
        c = colors{i};  m = markers{i};
        gain_curve = ami_opt_val(i,:) - ami_pam_val(i,:);
        plot(snr_vec, gain_curve, ['-' m], ...
            'Color', c, 'LineWidth', 2.0, ...
            'MarkerFaceColor', c, 'MarkerSize', 7, ...
            'DisplayName', sprintf('σ_R^2 = %.1f', sigX_vec(i)));
    end

    xlabel('SNR (dB)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Shaping Gain (bits/symbol)', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('GS Shaping Gain vs SNR — M=%d', M), ...
        'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 11);
    xticks(snr_vec);
    yline(0, 'k:', 'HandleVisibility', 'off');
    ylim([0, inf]);
    hold off;

    %% ====================================================================
    %  12a. FIGURE 3A: BER vs SNR (SEMILOG) — All channels
    % =====================================================================
    fig3a = figure('Name', 'BER vs SNR (Log Scale)', 'Color', 'w', ...
        'Units', 'normalized', 'Position', [0.02 0.05 0.45 0.40]);
    hold on; grid on; box on;

    for i = 1:nLevels
        c = colors{i};  m = markers{i};
        
        % Plot OPT (solid)
        semilogy(snr_vec, ber_opt(i,:), ['-' m], ...
            'Color', c, 'LineWidth', 2.0, ...
            'MarkerFaceColor', c, 'MarkerSize', 7, ...
            'DisplayName', sprintf('GS opt  (σ²_R=%.2f)', sigX_vec(i)));
        
        % Plot PAM (dashed)
        semilogy(snr_vec, ber_pam(i,:), ['--' m], ...
            'Color', c, 'LineWidth', 1.5, ...
            'MarkerFaceColor', 'w', 'MarkerSize', 6, ...
            'DisplayName', sprintf('PAM    (σ²_R=%.2f)', sigX_vec(i)));
    end

    xlabel('SNR (dB)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('BER', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('BER vs SNR — M=%d PAM (Log Scale, ML thresholds)', M), ...
        'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'southwest', 'NumColumns', 2, 'FontSize', 9);
    ylim([1e-10, 1]);
    xlim([snr_vec(1)-1, snr_vec(end)+1]);
    xticks(snr_vec);
    set(gca, 'YScale', 'log');
    hold off;

    %% ====================================================================
    %  12b. FIGURE 3B: BER vs SNR (LINEAR) — Turbulence Only
    % =====================================================================
    fig3b = figure('Name', 'BER vs SNR (Linear Scale)', 'Color', 'w', ...
        'Units', 'normalized', 'Position', [0.30 0.20 0.45 0.40]);
    hold on; grid on; box on;

    for i = 1:nLevels
        if sigX_vec(i) == 0, continue; end
        
        c = colors{i};  m = markers{i};
        
        % Plot OPT (solid)
        plot(snr_vec, ber_opt(i,:), ['-' m], ...
            'Color', c, 'LineWidth', 2.0, ...
            'MarkerFaceColor', c, 'MarkerSize', 7, ...
            'DisplayName', sprintf('GS opt  (σ²_R=%.2f)', sigX_vec(i)));
        
        % Plot PAM (dashed)
        plot(snr_vec, ber_pam(i,:), ['--' m], ...
            'Color', c, 'LineWidth', 1.5, ...
            'MarkerFaceColor', 'w', 'MarkerSize', 6, ...
            'DisplayName', sprintf('PAM    (σ²_R=%.2f)', sigX_vec(i)));
    end

    xlabel('SNR (dB)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('BER', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('BER vs SNR — M=%d PAM (Linear Scale, ML thresholds)', M), ...
        'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'northeast', 'NumColumns', 2, 'FontSize', 9);
    xlim([snr_vec(1)-1, snr_vec(end)+1]);
    xticks(snr_vec);
    hold off;

    %% ====================================================================
    %  13. FIGURE 4: SHAPING GAIN vs SNR (dB)
    % =====================================================================
    fig4 = figure('Name', 'Shaping Gain (dB) vs SNR', 'Color', 'w', ...
        'Units', 'normalized', 'Position', [0.50 0.05 0.45 0.40]);
    hold on; grid on; box on;

    for i = 1:nLevels
        if sigX_vec(i) == 0, continue; end
        c = colors{i};  m = markers{i};
        gain_dB_curve = 10*log10(ami_opt_val(i,:) ./ max(ami_pam_val(i,:), 1e-6));
        plot(snr_vec, gain_dB_curve, ['-' m], ...
            'Color', c, 'LineWidth', 2.0, ...
            'MarkerFaceColor', c, 'MarkerSize', 7, ...
            'DisplayName', sprintf('σ_R^2 = %.1f', sigX_vec(i)));
    end

    xlabel('SNR (dB)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Shaping Gain (dB)', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('GS Shaping Gain (dB) vs SNR — M=%d', M), ...
        'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 11);
    xticks(snr_vec);
    yline(0, 'k:', 'HandleVisibility', 'off');
    hold off;

    %% ====================================================================
    %  14. DEBUG LOG: OPTIMAL CONSTELLATIONS (for reproducibility)
    % =====================================================================
    fprintf('\n');
    fprintf('╔══════════════════════════════════════════════════════════════════════════════════════════════════════════════════╗\n');
    fprintf('║                              OPTIMAL CONSTELLATIONS (for reproducibility)                                        ║\n');
    fprintf('╚══════════════════════════════════════════════════════════════════════════════════════════════════════════════════╝\n\n');
    
    for i = 1:nLevels
        fprintf('── σ²_R = %.2f (%s) ──────────────────────────────────────────────────────────────────\n', ...
            sigX_vec(i), level_names(i));
        for j = 1:nSNR
            x_opt = x_opt_grid{i, j};
            fprintf('  SNR=%2ddB: x_opt = [', snr_vec(j));
            fprintf(' %.6f', x_opt);
            fprintf(' ]  (AMI=%.4f)\n', ami_opt_val(i,j));
        end
        fprintf('\n');
    end
    
    %% ====================================================================
    %  15. MATLAB CODE FOR RELOADING CONSTELLATIONS
    % =====================================================================
    fprintf('══════════════════════════════════════════════════════════════════════════════════════════════════════════════════\n');
    fprintf('MATLAB code to reproduce these constellations:\n');
    fprintf('══════════════════════════════════════════════════════════════════════════════════════════════════════════════════\n\n');
    
    fprintf('%% Optimal constellations found by SA (copy-paste to use)\n');
    fprintf('%% M=%d, P_avg=%.1f, power constraint = "mean" (IM/DD)\n\n', M, P_avg);
    
    for i = 1:nLevels
        fprintf('%% σ²_R = %.2f\n', sigX_vec(i));
        for j = 1:nSNR
            x_opt = x_opt_grid{i, j};
            fprintf('x_opt_M%d_sig%.0f_snr%02d = [', M, sigX_vec(i)*10, snr_vec(j));
            fprintf(' %.6f', x_opt);
            fprintf(' ]'';\n');
        end
        fprintf('\n');
    end
    
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

    cfg.sig_t = sqrt(log(1 + sigma_X_sq));
    cfg.mu_t  = -0.5 * cfg.sig_t^2;

    cfg.px = ones(M, 1) / M;

    x_max_bound = 5;
    cfg.ghN_h  = ghN_h;
    cfg.y_grid = AMI_functions.build_noCSI_y_grid(cfg, x_max_bound);

    sa = struct();
    sa.maxIter      = sa_maxIter;
    sa.itersPerTemp = 50;
    sa.T0           = 0.4;
    sa.Tf           = 1e-3;
    sa.nBlocks      = ceil(sa.maxIter / sa.itersPerTemp);
    sa.coolingRate  = exp(log(sa.Tf / sa.T0) / sa.nBlocks);

    sa.baseStd0      = 0.15;   sa.baseStdMin  = 1e-3;  sa.baseStdMax = 0.5;
    sa.targetAccLo   = 0.20;   sa.targetAccHi = 0.60;
    sa.baseStdGrow   = 1.25;   sa.baseStdShrink = 0.80;

    sa.enforce_sort    = true;
    sa.imdd_mode       = true;
    sa.powerConstraint = "mean";
    sa.enforce_power   = true;
    sa.minGap          = 0.05;
    sa.projectIters    = 2;

    sa.nStarts         = sa_nStarts;
    sa.useParallel     = sa_useParfor;
    sa.numWorkers      = [];
    sa.closePoolWhenDone = false;
    sa.logEvery        = 5000;

    cfg.SA = sa;

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