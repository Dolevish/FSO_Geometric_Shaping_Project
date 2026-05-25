function cfg = config_params_PhaseNoise()
% CONFIG_PARAMS_PHASENOISE
% Creates a configuration structure for the Phase Noise Optimization.
% 
% Defines:
%   1. System parameters (Constellation size, SNR, Phase Noise).
%   2. Simulated Annealing (SA) hyper-parameters.

    %% --- 1. System & Channel Parameters ---
    cfg.M        = 8;          % Modulation Order (e.g., 8 for 8-ary)
    cfg.P_avg    = 1;          % Average Power Constraint (E[|x|^2] <= P_avg)
    
    % Default Channel Conditions (Can be overridden in main script)
    cfg.SNR_dB   = 12;         % Signal-to-Noise Ratio [dB]
    cfg.PNSD_deg = 10;         % Phase Noise Standard Deviation [degrees]

    %% --- 2. Simulated Annealing (SA) Parameters ---
    % General SA Settings
    cfg.SA.maxIter      = 5000;   % Total number of iterations
    cfg.SA.printEvery   = 1000;    % Log progress every N iterations
    
    % Temperature Schedule (Cooling)
    cfg.SA.initialTemperature = 3; 
    cfg.SA.cooling_alpha      = 0.999; % T_new = T_old * alpha
    
    % Proposal Distribution (How we move points)
    cfg.SA.proposal_std_deg   = 10;    % Standard deviation for phase perturbation
    cfg.SA.amp_prob           = 0.1;   % Probability to choose Amplitude step over Phase step
    
    % Monte-Carlo Accuracy
    cfg.SA.Nmc_full     = 2e4;     % Number of MC samples for Mutual Information estimation
    cfg.SA.rngSeed_MI   = 1;       % Fixed seed for MI estimator (reduces noise in Cost Function)
    
    % Initialization
    % -1 = Random Start (Shuffle), >0 = Fixed Seed (Reproducible)
    cfg.SA.rngSeed_init = -1;      
end