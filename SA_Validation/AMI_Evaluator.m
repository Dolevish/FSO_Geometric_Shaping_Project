classdef AMI_Evaluator
    % AMI_EVALUATOR Evaluates Mutual Information for Phase Noise Channels.
    %
    % This class consolidates ALL MI calculations:
    %   1. Initialization (Monte-Carlo simulation)
    %   2. Incremental Updates (Fast Move - Log-Sum-Exp trick)
    %   3. Full Updates (Amplitude Change - consistency check)
    %
    % Usage:
    %   st = AMI_Evaluator(x0, cfg, Nmc, seed);  % Create State
    %   mi_new = st.eval_move(idx, u_new);       % Check Move
    %   st = st.apply_move(idx, u_new);          % Apply Move
    %   [mi, st] = st.recalc_full(x_new);        % Recalculate All
    
    properties
        % System Params
        M, N
        K_n, K_phi
        awgnOnly
        
        % Constellation & Channel State
        x           % Current Constellation (1xM)
        idxX        % Transmitted indices (1xN)
        phi         % Phase Noise realization (1xN)
        noise       % Thermal Noise realization (1xN)
        Y           % Received Signal (1xN)
        
        % Intermediate Calculations (for Incremental Updates)
        m           % Max Log-Likelihood per sample (1xN)
        S           % Sum of Exponentials per sample (1xN)
        logNum      % Numerator (Likelihood of sent symbol)
        logDen      % Denominator (Total Probability)
        
        % Precomputed constants
        log_const_awgn
        C
        logI0_Kphi
        penalty_u
        
        % Result
        mi_bits     % Current MI value
    end
    
    methods
        %% --- 1. CONSTRUCTOR (Initialization) ---
        function obj = AMI_Evaluator(x_init, params, Nmc, rngSeed)
            if nargin < 4 || isempty(rngSeed), rngSeed=1; end
            rng(rngSeed);
            
            % Setup Parameters
            obj.x = x_init(:).';
            obj.M = numel(obj.x);
            obj.N = Nmc;
            
            obj.K_n = 10^(params.SNR_dB / 10);
            sigma_n_sq = 1 / obj.K_n;
            
            pnsd_rad = params.PNSD_deg * pi / 180;
            if pnsd_rad < 1e-6
                obj.awgnOnly = true; obj.K_phi = Inf;
            else
                obj.awgnOnly = false; obj.K_phi = 1/(pnsd_rad^2);
            end
            
            % Monte-Carlo Simulation (Create the "World")
            obj.idxX = randi(obj.M, [1, Nmc]);
            
            if obj.awgnOnly
                obj.phi = zeros(1, Nmc);
            else
                if ~exist('circ_vmrnd','file')
                     obj.phi = (pnsd_rad * randn(1, Nmc)); 
                else
                     obj.phi = circ_vmrnd(0, obj.K_phi, Nmc).'; 
                end
            end
            
            obj.noise = (randn(1, Nmc) + 1i*randn(1, Nmc)) * sqrt(sigma_n_sq/2);
            
            % Calculate Y based on initial X
            obj = obj.update_channel_output();
            
            % Initial MI Calculation
            obj = obj.calc_mi_from_scratch();
        end
        
        %% --- 2. FAST INCREMENTAL UPDATE (Phase Steps) ---
        
        function mi_val = eval_move(obj, j, u_new)
            % EVAL_MOVE Fast "peek" at MI if we change symbol j -> u_new
            % Does NOT update the state. Returns the predicted MI.
            [mi_val, ~] = obj.core_update_logic(j, u_new, false);
        end
        
        function obj = apply_move(obj, j, u_new)
            % APPLY_MOVE Permanently apply change symbol j -> u_new
            % Updates all internal state matrices (m, S, etc.)
            [~, obj] = obj.core_update_logic(j, u_new, true);
        end
        
        %% --- 3. FULL RECALCULATION (Amplitude Steps) ---
        
        function [mi_val, obj] = recalc_full(obj, x_new)
            % RECALC_FULL Re-runs calculation on the SAME noise realization.
            % Used when global changes (like power normalization) occur.
            
            obj.x = x_new;
            obj = obj.update_channel_output(); % Update Y using stored noise
            obj = obj.calc_mi_from_scratch();  % Re-run full scanner
            mi_val = obj.mi_bits;
        end
        
    end
    
    methods (Access = private)
        
        %% --- Internal Logic Helpers ---
        
        function obj = update_channel_output(obj)
            % Recalculates Y based on current x, phi, and noise
            obj.Y = obj.x(obj.idxX) .* exp(1j*obj.phi) + obj.noise;
            
            % Update Precomputed Constants
            if obj.awgnOnly
                obj.log_const_awgn = log(obj.K_n/pi);
            else
                obj.C = log(obj.K_n/(2*pi)) - (obj.K_n/2)*abs(obj.Y).^2;
                obj.penalty_u = (-(obj.K_n/2) * abs(obj.x).^2).';
                obj.logI0_Kphi = obj.logI0_stable(obj.K_phi);
            end
        end
        
        function obj = calc_mi_from_scratch(obj)
            % The "Heavy" full calculator (like original AMI_Evaluator)
            m_curr = -inf(1, obj.N);
            S_curr = zeros(1, obj.N);
            logN = zeros(1, obj.N);
            
            for j=1:obj.M
                r_j = obj.row_loglik(obj.x(j));
                
                % Log-Sum-Exp Accumulation
                new_m = max(m_curr, r_j);
                S_curr = exp(m_curr - new_m).*S_curr + exp(r_j - new_m);
                m_curr = new_m;
                
                mask = (obj.idxX == j);
                if any(mask)
                    logN(mask) = r_j(mask);
                end
            end
            
            % Finalize
            obj.m = m_curr;
            obj.S = S_curr;
            obj.logNum = logN;
            obj.logDen = m_curr + log(S_curr);
            
            mi_nat = log(obj.M) + mean(obj.logNum - obj.logDen);
            obj.mi_bits = max(0, min(log2(obj.M), mi_nat/log(2)));
        end
        
        function [mi_new, obj_out] = core_update_logic(obj, j, u_new, do_apply)
            % The "Fast" Incremental Logic (moved from simulated_annealing)
            if do_apply, obj_out = obj; else, obj_out = obj; end
            
            m_loc = obj.m; 
            S_loc = obj.S; 
            logNum_loc = obj.logNum;
            
            mask_tx = (obj.idxX == j);
            
            % A. Update samples where X_j was transmitted
            if any(mask_tx)
                % Only need to update Y for these samples!
                Y_sub = u_new * exp(1j*obj.phi(mask_tx)) + obj.noise(mask_tx);
                
                % Calculate Likelihoods for ALL symbols vs these NEW Y values
                if obj.awgnOnly
                    r_all = obj.log_const_awgn - obj.K_n * abs(bsxfun(@minus, Y_sub, obj.x.')).^2;
                else
                    % Phase Noise math for subspace
                    C_sub = log(obj.K_n/(2*pi)) - (obj.K_n/2)*abs(Y_sub).^2;
                    pen   = obj.penalty_u; 
                    pen(j)= -(obj.K_n/2)*abs(u_new)^2; 
                    
                    Z = bsxfun(@times, conj(Y_sub), obj.x.');
                    alpha = obj.K_phi + obj.K_n*real(Z);
                    beta  = -obj.K_n*imag(Z);
                    arg   = sqrt(alpha.^2 + beta.^2);
                    
                    r_all = bsxfun(@plus, pen, (obj.logI0_stable(arg) - obj.logI0_Kphi)) + ...
                            repmat(C_sub, obj.M, 1);
                end
                
                m_loc(mask_tx) = max(r_all, [], 1);
                S_loc(mask_tx) = sum(exp(bsxfun(@minus, r_all, m_loc(mask_tx))), 1);
                logNum_loc(mask_tx) = r_all(j, :); % New numerator
                
                if do_apply
                    obj_out.Y(mask_tx) = Y_sub;
                    if ~obj.awgnOnly
                         % Update C and Penalty in state only if applied
                         % For strict correctness with cached C, we'd update obj.C(mask_tx) here.
                         obj_out.C(mask_tx) = log(obj.K_n/(2*pi)) - (obj.K_n/2)*abs(Y_sub).^2;
                         obj_out.penalty_u(j) = -(obj.K_n/2)*abs(u_new)^2;
                    end
                end
            end
            
            % B. Update samples where OTHER symbols were transmitted
            mask_nontx = ~mask_tx;
            if any(mask_nontx)
                u_old = obj.x(j);
                % Calc likelihood of OLD symbol j and NEW symbol j vs OLD Y
                r_old = obj.row_loglik_idx(u_old, mask_nontx);
                r_new = obj.row_loglik_idx(u_new, mask_nontx);
                
                mk = m_loc(mask_nontx);
                Sk = S_loc(mask_nontx);
                
                % Incremental Log-Sum-Exp Update
                inc = r_new > mk;
                same = ~inc;
                
                if any(same)
                    idx = find(same);
                    Sk(idx) = Sk(idx) - exp(r_old(idx)-mk(idx)) + exp(r_new(idx)-mk(idx));
                end
                if any(inc)
                    idx = find(inc);
                    Sk(idx) = exp(mk(idx)-r_new(idx)).*Sk(idx) - exp(r_old(idx)-r_new(idx)) + 1;
                    mk(idx) = r_new(idx);
                end
                
                m_loc(mask_nontx) = mk;
                S_loc(mask_nontx) = Sk;
            end
            
            % C. Final MI
            logDen_loc = m_loc + log(S_loc);
            mi_nat = log(obj.M) + mean(logNum_loc - logDen_loc);
            mi_new = max(0, min(log2(obj.M), mi_nat/log(2)));
            
            if do_apply
                obj_out.m = m_loc; obj_out.S = S_loc; 
                obj_out.logNum = logNum_loc; obj_out.logDen = logDen_loc;
                obj_out.mi_bits = mi_new;
                obj_out.x(j) = u_new;
            end
        end
        
        function r = row_loglik(obj, u)
            % Calculates log P(Y|u) for ALL samples
            if obj.awgnOnly
                r = obj.log_const_awgn - obj.K_n * (abs(obj.Y - u).^2);
            else
                z = conj(obj.Y) .* u;
                alpha = obj.K_phi + obj.K_n * real(z);
                beta  = -obj.K_n * imag(z);
                arg   = sqrt(alpha.^2 + beta.^2);
                penalty = -(obj.K_n/2)*abs(u)^2;
                r = obj.C - obj.logI0_Kphi + penalty + obj.logI0_stable(arg);
            end
        end
        
        function r = row_loglik_idx(obj, u, mask)
            % Calculates log P(Y|u) for SUBSET of samples
            Y_sub = obj.Y(mask);
            if obj.awgnOnly
                r = obj.log_const_awgn - obj.K_n * (abs(Y_sub - u).^2);
            else
                z = conj(Y_sub) .* u;
                alpha = obj.K_phi + obj.K_n * real(z);
                beta  = -obj.K_n * imag(z);
                arg   = sqrt(alpha.^2 + beta.^2);
                penalty = -(obj.K_n/2)*abs(u)^2;
                r = obj.C(mask) - obj.logI0_Kphi + penalty + obj.logI0_stable(arg);
            end
        end

    end
    
    methods (Static)
        function v = logI0_stable(x)
            v = log(max(besseli(0, x, 1), realmin)) + x;
        end
    end
end