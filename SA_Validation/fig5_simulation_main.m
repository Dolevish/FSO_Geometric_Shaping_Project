function fig5_simulation_main



    clear; close all; clc;

    rng('shuffle');

    % --- experiment settings ---
    snr_vec  = 1:2:13;
    pnsd_vec = [0, 10, 25, 40];
    num_random_retries = 4;
    
    % flatten (SNR, PNSD) grid for parfor
    [SNR_Grid, PNSD_Grid] = meshgrid(snr_vec, pnsd_vec);
    task_snr  = SNR_Grid(:);
    task_pnsd = PNSD_Grid(:);
    num_tasks = numel(task_snr);
    
    total_optimizations = num_tasks * (1 + num_random_retries);
    
    % --- base configuration ---
    base_cfg = config_params_PhaseNoise();
    base_cfg.M = 8;
    base_cfg.P_avg = 1;
    
    % SA settings
    base_cfg.SA.maxIter    = 3000;
    base_cfg.SA.Nmc_full   = 2e4;
    base_cfg.SA.printEvery = 1500;
    
    % --- constellations ---
    x_8psk = define_constellation_2D(base_cfg.M, base_cfg.P_avg, 'PSK');
    
    % --- linear result buffers (for parfor) ---
    results_opt_linear = zeros(num_tasks, 1);
    results_psk_linear = zeros(num_tasks, 1);

    % pre-generate random seeds for MI Monte-Carlo
    mi_seeds = randi(2^31-1, num_tasks, 1);
    
    % --- parallel pool creation ---
    pool = gcp('nocreate');
    if isempty(pool)
        pool = parpool('local');
    end
    numWorkers = pool.NumWorkers;
    
    % --- progress counter ---
    Q = parallel.pool.DataQueue;
    afterEach(Q, @update_progress);
    fprintf('=== Starting SA Validation (%d tasks on %d workers) ===\n', ...
            num_tasks, numWorkers);
    
    total_start = tic;
    
	parfor k = 1:num_tasks
		current_snr  = task_snr(k);
		current_pnsd = task_pnsd(k);
		
		% Create unique label for this worker task
		taskLabel = sprintf('[Task=%d|PNSD=%d|SNR=%d]', k, current_pnsd,current_snr);
		
		cfg = base_cfg;
		cfg.PNSD_deg = current_pnsd;
		cfg.SNR_dB   = current_snr;

		% 8-PSK reference with random seed per task
		st_ref = AMI_Evaluator(x_8psk, cfg, 5e4, mi_seeds(k));
		results_psk_linear(k) = st_ref.mi_bits;
		
		% SA optimization
		best_mi_point = -inf;
		
		% QAM-based init
		cfg.SA.rngSeed_init = -1; 
		try
			[~, mi_qam, ~] = simulated_annealing(cfg, x_8psk, taskLabel);
		catch
			[~, mi_qam, ~] = simulated_annealing(cfg, [], taskLabel); 
		end
		if mi_qam > best_mi_point, best_mi_point = mi_qam; end
		send(Q, 1);
		
		% Random restarts
		for r = 1:num_random_retries

			[~, mi_rnd, ~] = simulated_annealing(cfg, [], taskLabel);
			if mi_rnd > best_mi_point, best_mi_point = mi_rnd; end
			send(Q, 1);
		end
		
		results_opt_linear(k) = best_mi_point;
	end
	
    total_time = toc(total_start);
    fprintf('\nDone! Total time: %.1f minutes.\n', total_time/60);
    
    ami_opt = reshape(results_opt_linear, size(SNR_Grid));
    ami_psk = reshape(results_psk_linear, size(SNR_Grid));
    
    % plot results
    plot_final_results(snr_vec, pnsd_vec, ami_opt, ami_psk);
end

function update_progress(~)
    % parallel progress counter
    persistent count
    if isempty(count), count = 0; end
    count = count + 1;
    if mod(count, 5) == 0
        fprintf('Progress: %d runs completed.\n', count);
    end
end

function plot_final_results(snr_vec, pnsd_vec, ami_opt, ami_psk)
    figure('Name', 'Final Hybrid Results (Optimized)', 'Color', 'w', 'Position', [100 100 900 700]);
    hold on; grid on; box on;
    
    colors = {
        [0 0.4470 0.7410], ... 
        [0.6350 0.0780 0.1840], ... 
        [0.4660 0.6740 0.1880], ... 
        [0.9290 0.6940 0.1250]      
    };
    
    markers = {'o', 's', '^', 'd'};
    
    for i = 1:length(pnsd_vec)
        c = colors{i};
        m = markers{i};
        pnsd = pnsd_vec(i);
        
        plot(snr_vec, ami_opt(i, :), ['-' m], 'Color', c, 'LineWidth', 2, ...
             'MarkerFaceColor', c, 'MarkerSize', 7, ...
             'DisplayName', sprintf('Optimized (PNSD=%d)', pnsd));
         
        plot(snr_vec, ami_psk(i, :), ['--' m], 'Color', c, 'LineWidth', 1.5, ...
             'MarkerFaceColor', 'w', 'MarkerSize', 6, ...
             'DisplayName', sprintf('8-PSK (PNSD=%d)', pnsd));
    end
    
    yline(3, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off'); 
    
    xlabel('SNR (dB)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('AMI (bits/symbol)', 'FontSize', 12, 'FontWeight', 'bold');
    title('AMI Performance: SA Optimized vs 8-PSK', 'FontSize', 14);
    
    legend('Location', 'SouthEast', 'NumColumns', 2, 'FontSize', 10);
    ylim([0.5, 3.1]);
    xlim([1, 13]);
    hold off;
end
