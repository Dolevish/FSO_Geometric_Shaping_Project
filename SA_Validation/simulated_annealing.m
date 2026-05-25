function [x_best, mi_best, info] = simulated_annealing(cfg, x_init, taskLabel)
% Algorithm: Simulated Annealing .
%
% INPUTS:
%   cfg    - Configuration struct (from config_params_PhaseNoise)
%   x_init - (Optional) Initial constellation vector. If empty, random init.
%
% OUTPUTS:
%   x_best  - Optimized constellation points
%   mi_best - The maximum Mutual Information achieved
%   info    - Struct with convergence history and timing

    %% --- 1. Initialization ---
    if nargin < 1 || isempty(cfg)
        cfg = config_params_PhaseNoise();
    end
    
    sa = cfg.SA;
    M  = cfg.M;
    P_avg = cfg.P_avg;

    % Handle missing arguments
    if nargin < 3, taskLabel = '[SA]'; end
    if nargin < 2, x_init = []; end

    if ~isempty(x_init)
        if numel(x_init) ~= M
            error('simulated_annealing: x_init must have exactly M points.');
        end
        x0 = x_init(:);
        fprintf('%s Initializing with CUSTOM constellation.\n', taskLabel);
        rng('shuffle'); 
    else
        % Random Initialization
        if isfield(sa, 'rngSeed_init') && sa.rngSeed_init > 0
            rng(sa.rngSeed_init);
            fprintf('%s Initializing with FIXED seed: %d\n', taskLabel, sa.rngSeed_init);
        else
            rng('shuffle');
            fprintf('%s Initializing with RANDOM seed...\n', taskLabel);
        end
        x0 = (randn(M,1) + 1i*randn(M,1));
    end
    
    % Enforce power constraint initially
    x0 = Enforce_Power_Constraint(x0, P_avg);

    % Initial Cost Evaluation (Maximize MI)
    st = AMI_Evaluator(x0, cfg, sa.Nmc_full, sa.rngSeed_MI);
    
    x_cur  = st.x;
    mi_cur = st.mi_bits;
    
    x_best  = x_cur;
    mi_best = mi_cur;
    
    T = sa.initialTemperature;
    
    % History logging
    info.mi_hist = nan(sa.maxIter, 1);
    info.T_hist  = nan(sa.maxIter, 1);
    
    t_start = tic;

    %% --- 2. Main Optimization Loop ---
    for it = 1:sa.maxIter
        
        % Propose a Move
        j = randi(M); % Pick a random symbol to modify
        
        % Decide move type: Amplitude change or Phase rotation
        do_amplitude_step = (rand() < sa.amp_prob); 
        
        if do_amplitude_step
            % -- Amplitude Step (Global Normalization) --
            r_scale = 1 + 0.05 * randn(); 
            x_new_temp = x_cur;
            x_new_temp(j) = x_new_temp(j) * r_scale;
            
            % Global normalization affects ALL points
            x_new = Enforce_Power_Constraint(x_new_temp, P_avg);
            
            [mi_new, st_new_amp] = st.recalc_full(x_new);
            
        else
            % -- Phase Step--
            dphi = (sa.proposal_std_deg * pi/180) * randn();
            u_prop = x_cur(j) * exp(1j * dphi);
            
            mi_new = st.eval_move(j, u_prop);
            
            % Virtual update variables
            x_new = x_cur; 
            x_new(j) = u_prop;     
        end

        dE = -(mi_new - mi_cur); 

        accept = false;
        if dE <= 0
            accept = true; % Always accept improvement
        else
            T_eff = max(T, 1e-9);
            prob = exp(-dE / T_eff);
            if rand() < prob
                accept = true;
            end
        end

        % Update State
        if accept
            mi_cur = mi_new;
            
            if do_amplitude_step
                x_cur = x_new;
                st = st_new_amp; 
            else

                st = st.apply_move(j, x_new(j));
                x_cur = st.x;
            end

            % Track Best Solution
            if mi_cur > mi_best
                mi_best = mi_cur;
                x_best  = x_cur;
            end
        end

        % D. Cooling & Logging
        T = T * sa.cooling_alpha;
        info.mi_hist(it) = mi_cur;
        info.T_hist(it)  = T;

        % Print progress
        if it == 1 || mod(it, sa.printEvery) == 0 || it == sa.maxIter
             elapsed = toc(t_start);
             % This function is defined below as a local function
             eta_str = get_eta_str(it, sa.maxIter, elapsed);
             
             fprintf('%s [SA] Iter %5d/%d | T=%.1e | MI=%.5f | Best=%.5f | ETA: %s\n', ...
                     taskLabel, it, sa.maxIter, T, mi_cur, mi_best, eta_str);
             
             % debug print
            fprintf('    x_cur  = %s\n', fmt_vec_complex(x_cur));
            fprintf('    x_best = %s\n', fmt_vec_complex(x_best));
        end
    end

    info.time = toc(t_start);
    info.x0 = x0;
end

%% --- Local Helper Functions ---

function s = get_eta_str(it, maxIter, elapsed)
% GET_ETA_STR Calculates estimated time remaining
    if it <= 1 || elapsed <= 0
        s = '...';
    else
        rem_sec = (maxIter - it) * (elapsed / it);
        if rem_sec < 60
            s = sprintf('%.0fs', rem_sec);
        else
            s = sprintf('%.1fm', rem_sec/60);
        end
    end
end

function s = fmt_vec_complex(v)
% FMT_VEC_COMPLEX Formats complex vector for printing
    s = sprintf('[%s]', strjoin(arrayfun(@(z) sprintf('%+.2f%+.2fi', real(z), imag(z)), v(:).', 'UniformOutput', false), ', '));
end