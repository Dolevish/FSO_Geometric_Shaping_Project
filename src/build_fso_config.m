function cfg = build_fso_config(M, P_avg, SNR_dB, sigma_X_sq, varargin)
%BUILD_FSO_CONFIG  Canonical configuration for the No-CSI FSO GS project.
%
%   cfg = build_fso_config(M, P_avg, SNR_dB, sigma_X_sq)
%   cfg = build_fso_config(..., 'Name', Value, ...)
%
% Core channel model:
%       y = h R x + n
% with h log-normal, E[h]=1, and n~N(0,sigma_n^2).
%
% The SNR convention follows the manuscript definition
%       SNR = P_avg^2 / sigma_n^2
% so
%       sigma_n^2 = P_avg^2 / 10^(SNR_dB/10).
%
% sigma_X_sq is the normalized intensity variance / scintillation-index
% parameter used by the code.  With E[h]=1,
%       sigma_X_sq = Var(h).
%
% Name-value options:
%   'R'              detector responsivity (default 1)
%   'ghN_h'          GH order for fast fading integral (default 40)
%   'xMaxBound'      bound used to pre-build fast y-grid (default 5)
%   'saMaxIter'      SA iterations per start (default 10000)
%   'saNStarts'      number of SA restarts (default 8)
%   'saUseParallel'  enable parallel restarts (default false)
%   'minGap'         adjacent-level minimum spacing (default 0.05)
%   'pinZero'        enforce x(1)=0 in IM/DD projection (default true)
%   'seedInit'       master RNG seed; [] means random (default [])
%   'logEvery'       SA progress-print interval (default 5000)
%
% This file is the single source of truth for channel and SA defaults in the
% IEEE revision.  Simulation drivers should migrate to this builder instead
% of maintaining local build_config() copies.

    p = inputParser;
    p.FunctionName = mfilename;

    addRequired(p, 'M', @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v >= 2 && mod(v,1)==0);
    addRequired(p, 'P_avg', @(v) isnumeric(v) && isscalar(v) && isfinite(v) && isreal(v) && v > 0);
    addRequired(p, 'SNR_dB', @(v) isnumeric(v) && isscalar(v) && isfinite(v) && isreal(v));
    addRequired(p, 'sigma_X_sq', @(v) isnumeric(v) && isscalar(v) && isfinite(v) && isreal(v) && v >= 0);

    addParameter(p, 'R', 1, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && isreal(v) && v > 0);
    addParameter(p, 'ghN_h', 40, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v >= 2 && mod(v,1)==0);
    addParameter(p, 'xMaxBound', 5, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && isreal(v) && v > 0);

    addParameter(p, 'saMaxIter', 10000, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v >= 1 && mod(v,1)==0);
    addParameter(p, 'saNStarts', 8, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v >= 1 && mod(v,1)==0);
    addParameter(p, 'saUseParallel', false, @(v) islogical(v) && isscalar(v));
    addParameter(p, 'minGap', 0.05, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && isreal(v) && v >= 0);
    addParameter(p, 'pinZero', true, @(v) islogical(v) && isscalar(v));
    addParameter(p, 'seedInit', [], @(v) isempty(v) || (isnumeric(v) && isscalar(v) && isfinite(v) && v >= 0));
    addParameter(p, 'logEvery', 5000, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v >= 0 && mod(v,1)==0);

    parse(p, M, P_avg, SNR_dB, sigma_X_sq, varargin{:});
    o = p.Results;

    cfg = struct();

    % ---------------------------------------------------------------------
    % Channel / constellation parameters
    % ---------------------------------------------------------------------
    cfg.M          = double(M);
    cfg.P_avg      = double(P_avg);
    cfg.SNR_dB     = double(SNR_dB);
    cfg.sigma_X_sq = double(sigma_X_sq);
    cfg.R          = double(o.R);

    SNR_lin        = 10^(cfg.SNR_dB / 10);
    cfg.sigma_n_sq = cfg.P_avg^2 / SNR_lin;

    % Log-normal fading parametrization with E[h]=1.
    sig_t_sq = log(1 + cfg.sigma_X_sq);
    cfg.sig_t = sqrt(sig_t_sq);
    cfg.mu_t  = -0.5 * sig_t_sq;

    cfg.px = ones(cfg.M, 1) / cfg.M;

    % ---------------------------------------------------------------------
    % Fast AMI evaluator settings
    % ---------------------------------------------------------------------
    cfg.ghN_h     = double(o.ghN_h);
    cfg.xMaxBound = double(o.xMaxBound);
    cfg.y_grid    = AMI_functions.build_noCSI_y_grid(cfg, cfg.xMaxBound);

    % ---------------------------------------------------------------------
    % Simulated Annealing settings
    % ---------------------------------------------------------------------
    sa = struct();

    sa.maxIter      = double(o.saMaxIter);
    sa.itersPerTemp = 50;
    sa.T0           = 0.4;
    sa.Tf           = 1e-3;
    sa.nBlocks      = ceil(sa.maxIter / sa.itersPerTemp);
    sa.coolingRate  = exp(log(sa.Tf / sa.T0) / sa.nBlocks);

    sa.baseStd0      = 0.15;
    sa.baseStdMin    = 1e-3;
    sa.baseStdMax    = 0.5;
    sa.targetAccLo   = 0.20;
    sa.targetAccHi   = 0.60;
    sa.baseStdGrow   = 1.25;
    sa.baseStdShrink = 0.80;

    sa.enforce_sort    = true;
    sa.imdd_mode       = true;
    sa.pinZero         = logical(o.pinZero);
    sa.powerConstraint = "mean";
    sa.enforce_power   = true;
    sa.minGap          = double(o.minGap);

    % Kept only for backward compatibility with older scripts.  The new
    % projection is deterministic/idempotent and does not iterate projection.
    sa.projectIters = 1;

    sa.nStarts           = double(o.saNStarts);
    sa.useParallel       = logical(o.saUseParallel);
    sa.numWorkers        = [];
    sa.closePoolWhenDone = false;
    sa.logEvery          = double(o.logEvery);
    sa.seedInit          = o.seedInit;

    cfg.SA = sa;

    % Fail early if the requested min-gap is impossible under mean power.
    minRequiredMean = (cfg.M - 1) * cfg.SA.minGap / 2;
    if cfg.P_avg < minRequiredMean - 100*eps(max(1,cfg.P_avg))
        error('build_fso_config:InfeasibleMinGap', ...
            ['M=%d and minGap=%.6g require mean optical power at least %.6g, ' ...
             'but P_avg=%.6g.'], ...
            cfg.M, cfg.SA.minGap, minRequiredMean, cfg.P_avg);
    end

    % IMPORTANT: the evaluator is intentionally side-effect free.  Candidate
    % projection is the optimizer's responsibility; the evaluator only scores
    % the constellation it receives.
    params = cfg;
    y_grid_local = cfg.y_grid;
    ghN_local    = cfg.ghN_h;
    cfg.AMI_Evaluator = @(x_in) AMI_functions.AMI_noCSI_fast_grid( ...
        x_in, params.px, params, ghN_local, y_grid_local);
end
