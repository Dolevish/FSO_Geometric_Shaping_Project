function fig2_simulation_main

    clear;
    close all;
    clc;
    addpath(pwd);

    fprintf('Starting SA Validation (Kayhan & Montorsi, Fig 2) ...\n');

    % number of random SA initializations per PNSD
    num_random_inits = 4;

    rng('shuffle');

    base = config_params_PhaseNoise();
    
    % PNSD values
    PNSD_list = [0 10 25 40];
    nP = numel(PNSD_list);

    % result vectors
    mi_SA       = zeros(1, nP);
    mi_8psk     = zeros(1, nP); 
    mi_rand_vec = zeros(1, nP); 
    
    x_SA   = cell(1, nP);
    x0_SA  = cell(1, nP);

    % 8-PSK reference constellation
    x_8psk_const = define_constellation_2D(base.M, base.P_avg, 'PSK');
    
    Nmc_ref = 5e4;

    % --- 8-PSK reference---
    for k = 1:nP
        params_k = base;
        params_k.PNSD_deg = PNSD_list(k);
        fprintf('Calculating 8-PSK Reference for PNSD=%d...\n', params_k.PNSD_deg);
        st_ref = AMI_Evaluator(x_8psk_const, params_k, Nmc_ref, 123);
        mi_8psk(k) = st_ref.mi_bits;
    end

    % --- build task grid---
    [P_idx_grid, R_idx_grid] = ndgrid(1:nP, 1:num_random_inits);
    task_P_idx = P_idx_grid(:);      % which PNSD for each task
    task_R_idx = R_idx_grid(:);      % which random start index 
    num_tasks  = numel(task_P_idx);

    % linear buffers for all tasks
    mi_opt_lin  = zeros(num_tasks, 1);
    mi_rand_lin = zeros(num_tasks, 1);
    x_opt_lin   = cell(num_tasks, 1);
    x0_lin      = cell(num_tasks, 1);

    % random seeds for MI of random starts
    mi_rand_seeds = randi(2^31-1, num_tasks, 1);

    % --- parallel pool---
    pool = gcp('nocreate');
    if isempty(pool)
        pool = parpool('local');
    end
    numWorkers = pool.NumWorkers;
    fprintf('Using %d workers.\n', numWorkers);

    fprintf('Running %d tasks (%d PNSD values x %d random starts) in parallel...\n', ...
            num_tasks, nP, num_random_inits);


    parfor t = 1:num_tasks
        k = task_P_idx(t);  % PNSD index

        % local params
        params_local = base;
        params_local.PNSD_deg = PNSD_list(k);

        % Create unique label for this worker task
        taskLabel = sprintf('[Task=%d|PNSD=%d]', t, PNSD_list(k));

        % SA optimization from random start (pass [] for x_init)
        [x_opt_t, mi_opt_t, info_t] = simulated_annealing(params_local, [], taskLabel);

        % MI of the random-start constellation
        seed_t   = mi_rand_seeds(t);
        st_rand_t = AMI_Evaluator(info_t.x0, params_local, Nmc_ref, seed_t);
        mi_rand_t = st_rand_t.mi_bits;

        mi_opt_lin(t)  = mi_opt_t;
        mi_rand_lin(t) = mi_rand_t;
        x_opt_lin{t}   = x_opt_t;
        x0_lin{t}      = info_t.x0;
    end

    % --- pick best run per PNSD ---
    for k = 1:nP
        task_mask = (task_P_idx == k);
        idx_tasks = find(task_mask);

        mi_opt_k  = mi_opt_lin(idx_tasks);
        mi_rand_k = mi_rand_lin(idx_tasks);

        [mi_best, idx_local_best] = max(mi_opt_k);
        idx_best = idx_tasks(idx_local_best);

        mi_SA(k)       = mi_best;
        mi_rand_vec(k) = mi_rand_k(idx_local_best);
        x_SA{k}        = x_opt_lin{idx_best};
        x0_SA{k}       = x0_lin{idx_best};

        gain = mi_SA(k) - mi_8psk(k);
        fprintf('>>> BEST RESULT: PNSD=%d -> SA=%.4f, Rand=%.4f, 8PSK=%.4f (Gain=%.3f)\n', ...
                PNSD_list(k), mi_SA(k), mi_rand_vec(k), mi_8psk(k), gain);

        % best constellation plot for this PNSD
        figure('Name', sprintf('Optimization Result PNSD %d (best of %d)', ...
                PNSD_list(k), num_random_inits));
        plot(real(x0_SA{k}), imag(x0_SA{k}), 'bo', 'MarkerSize', 6, 'LineWidth', 1.0); hold on;
        plot(real(x_SA{k}),  imag(x_SA{k}),  'rx', 'MarkerSize', 9, 'LineWidth', 2.0);
        
        theta = linspace(0, 2*pi, 100);
        plot(cos(theta)*sqrt(base.P_avg), sin(theta)*sqrt(base.P_avg), 'k:', 'LineWidth', 0.5);
        
        axis equal; grid on; box on;
        xlabel('In-phase'); ylabel('Quadrature');
        title(sprintf('PNSD=%d | SA=%.3f | Rand=%.3f | 8PSK=%.3f (best of %d)', ...
              PNSD_list(k), mi_SA(k), mi_rand_vec(k), mi_8psk(k), num_random_inits));
        legend('Start (Rand)','Optimized','Unit Circle','Location','bestoutside');
        
        drawnow;
    end

    % --- summary figure ---
    figure('Name', 'MI Summary');
    
    plot(PNSD_list, mi_8psk, 'ks-', 'LineWidth', 1.5, 'MarkerSize', 8, ...
         'DisplayName', '8-PSK (Ref)');
    hold on;
    
    plot(PNSD_list, mi_rand_vec, 'b^--', 'LineWidth', 1.5, 'MarkerSize', 8, ...
         'DisplayName', 'Random (Start, best)');
    
    plot(PNSD_list, mi_SA, 'ro-', 'LineWidth', 2.0, 'MarkerSize', 8, ...
         'DisplayName', 'SA Optimized (best)');
    
    grid on; box on;
    xlabel('PNSD [deg]'); ylabel('MI [bits/sym]');
    title(sprintf('MI vs Phase Noise (SNR=%.1f dB) [ %d starts]', ...
          base.SNR_dB, num_random_inits));
    legend('Location','SouthWest');
    
    text(PNSD_list(end), mi_SA(end), ...
         sprintf('  Gain: +%.2f', mi_SA(end)-mi_8psk(end)), ...
         'Color', 'r', 'FontWeight', 'bold');
end
