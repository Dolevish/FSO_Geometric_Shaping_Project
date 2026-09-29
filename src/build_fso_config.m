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
% parameter used by the code. With E[h]=1,
%       sigma_X_sq = Var(h).
%
% Name-value options:
%   'R'                 detector responsivity (default 1)
%   'FastFadingMethod'  fast fading quadrature: 'logh' or 'gh'
%                       (default 'logh' from Commit I.3)
%   'ghN_h'             GH order when FastFadingMethod='gh' (default 40)
%   'loghDt'            uniform t=ln(h) spacing (default 0.01)
%   'loghSpanSigma'     integration half-span in sigma_t (default 8)
%   'loghYBlockSize'    y-block size for bounded memory (default 512)
%   'YGridMode'         'fixed-bound' or 'candidate-adaptive'
%                       (default 'fixed-bound'; Commit I.6 keeps historical
%                       behavior until optimizer-level A/B validation passes)
%   'YGridScale'        adaptive support multiplier >=1 (default 1.1)
%   'xMaxBound'         bound used to pre-build fixed y-grid (default 5)
%   'saMaxIter'         SA iterations per start (default 10000)
%   'saNStarts'         number of SA restarts (default 8)
%   'saUseParallel'     enable parallel restarts (default false)
%   'minGap'            adjacent-level minimum spacing (default 0.05)
%   'pinZero'           enforce x(1)=0 in IM/DD projection (default true)
%   'seedInit'          master RNG seed; [] means random (default [])
%   'logEvery'          SA progress-print interval (default 5000)
%   'historyEvery'      sample SA convergence history every N iterations;
%                       [] means one sample per temperature block (default [])
%
% The canonical configuration exposes two side-effect-free evaluators:
%   cfg.AMI_Evaluator  - fast objective used inside SA.
%   cfg.AMI_Validator  - independent validation evaluator used after SA.
%
% Commit I.6 y-grid semantics:
%   fixed-bound         : build y-grid once from cfg.xMaxBound and reuse it.
%   candidate-adaptive  : for each candidate x, rebuild y-grid from
%                         cfg.yGridScale*max(abs(x)). This changes only the
%                         numerical integration support, never the feasible
%                         constellation set and never introduces clipping.
%
% candidate-adaptive is currently supported only with the production log-h
% fading evaluator. GH remains fixed-bound to preserve Commit-G/I.4 meaning.
%
% To reproduce the pre-I.3 fast objective explicitly use:
%       'FastFadingMethod','gh','ghN_h',<order>
%
% This file is the single source of truth for channel and SA defaults in the
% IEEE revision. Simulation drivers should migrate to this builder instead
% of maintaining local build_config() copies.

    p = inputParser;
    p.FunctionName = mfilename;

    addRequired(p, 'M', @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v >= 2 && mod(v,1)==0);
    addRequired(p, 'P_avg', @(v) isnumeric(v) && isscalar(v) && isfinite(v) && isreal(v) && v > 0);
    addRequired(p, 'SNR_dB', @(v) isnumeric(v) && isscalar(v) && isfinite(v) && isreal(v));
    addRequired(p, 'sigma_X_sq', @(v) isnumeric(v) && isscalar(v) && isfinite(v) && isreal(v) && v >= 0);

    addParameter(p, 'R', 1, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && isreal(v) && v > 0);
    addParameter(p, 'FastFadingMethod', 'logh', @(v) ischar(v) || (isstring(v) && isscalar(v)));
    addParameter(p, 'ghN_h', 40, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v >= 2 && mod(v,1)==0);
    addParameter(p, 'loghDt', 0.01, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && isreal(v) && v > 0);
    addParameter(p, 'loghSpanSigma', 8, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && isreal(v) && v > 0);
    addParameter(p, 'loghYBlockSize', 512, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v >= 1 && mod(v,1)==0);
    addParameter(p, 'YGridMode', 'fixed-bound', @(v) ischar(v) || (isstring(v) && isscalar(v)));
    addParameter(p, 'YGridScale', 1.1, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && isreal(v) && v >= 1);
    addParameter(p, 'xMaxBound', 5, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && isreal(v) && v > 0);

    addParameter(p, 'saMaxIter', 10000, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v >= 1 && mod(v,1)==0);
    addParameter(p, 'saNStarts', 8, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v >= 1 && mod(v,1)==0);
    addParameter(p, 'saUseParallel', false, @(v) islogical(v) && isscalar(v));
    addParameter(p, 'minGap', 0.05, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && isreal(v) && v >= 0);
    addParameter(p, 'pinZero', true, @(v) islogical(v) && isscalar(v));
    addParameter(p, 'seedInit', [], @(v) isempty(v) || (isnumeric(v) && isscalar(v) && isfinite(v) && v >= 0));
    addParameter(p, 'logEvery', 5000, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v >= 0 && mod(v,1)==0);
    addParameter(p, 'historyEvery', [], @(v) isempty(v) || (isnumeric(v) && isscalar(v) && isfinite(v) && v >= 1 && mod(v,1)==0));

    parse(p, M, P_avg, SNR_dB, sigma_X_sq, varargin{:});
    o = p.Results;

    fastMethod = lower(string(o.FastFadingMethod));
    if ~ismember(fastMethod, ["logh","gh"])
        error('build_fso_config:BadFastFadingMethod', ...
            'FastFadingMethod must be ''logh'' or ''gh'', got ''%s''.', fastMethod);
    end

    yGridMode = lower(string(o.YGridMode));
    if ~ismember(yGridMode, ["fixed-bound","candidate-adaptive"])
        error('build_fso_config:BadYGridMode', ...
            'YGridMode must be ''fixed-bound'' or ''candidate-adaptive'', got ''%s''.', yGridMode);
    end
    if fastMethod == "gh" && yGridMode == "candidate-adaptive"
        error('build_fso_config:AdaptiveYGridRequiresLogh', ...
            ['candidate-adaptive y-grid is intentionally enabled only for the log-h ' ...
             'fast evaluator. Use YGridMode=''fixed-bound'' with GH.']);
    end

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
    cfg.fastFadingMethod = char(fastMethod);
    cfg.ghN_h            = double(o.ghN_h);
    cfg.loghDt           = double(o.loghDt);
    cfg.loghSpanSigma    = double(o.loghSpanSigma);
    cfg.loghYBlockSize   = double(o.loghYBlockSize);
    cfg.yGridMode        = char(yGridMode);
    cfg.yGridScale       = double(o.YGridScale);
    cfg.xMaxBound        = double(o.xMaxBound);

    % Always retain the historical fixed-bound grid in cfg for diagnostics,
    % reproducibility, and config snapshots. In candidate-adaptive mode the
    % AMI evaluator does not use this grid; it rebuilds support from x.
    cfg.y_grid = AMI_functions.build_noCSI_y_grid(cfg, cfg.xMaxBound);

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

    % Kept only for backward compatibility with older scripts. The new
    % projection is deterministic/idempotent and does not iterate projection.
    sa.projectIters = 1;

    sa.nStarts           = double(o.saNStarts);
    sa.useParallel       = logical(o.saUseParallel);
    sa.numWorkers        = [];
    sa.closePoolWhenDone = false;
    sa.logEvery          = double(o.logEvery);
    sa.seedInit          = o.seedInit;

    % Commit B: convergence/observability sampling cadence. Empty means
    % simulated_annealing samples once per temperature block.
    if isempty(o.historyEvery)
        sa.historyEvery = sa.itersPerTemp;
    else
        sa.historyEvery = double(o.historyEvery);
    end

    cfg.SA = sa;

    % Fail early if the requested min-gap is impossible under mean power.
    minRequiredMean = (cfg.M - 1) * cfg.SA.minGap / 2;
    if cfg.P_avg < minRequiredMean - 100*eps(max(1,cfg.P_avg))
        error('build_fso_config:InfeasibleMinGap', ...
            ['M=%d and minGap=%.6g require mean optical power at least %.6g, ' ...
             'but P_avg=%.6g.'], ...
            cfg.M, cfg.SA.minGap, minRequiredMean, cfg.P_avg);
    end

    % IMPORTANT: both evaluators are side-effect free. Candidate projection
    % is the optimizer's responsibility; evaluators only score the exact
    % constellation they receive.
    params = cfg;
    y_grid_local = cfg.y_grid;

    switch fastMethod
        case "gh"
            ghN_local = cfg.ghN_h;
            cfg.AMI_Evaluator = @(x_in) AMI_functions.AMI_noCSI_fast_grid( ...
                x_in, params.px, params, ghN_local, y_grid_local);

        case "logh"
            dt_local = cfg.loghDt;
            span_local = cfg.loghSpanSigma;
            block_local = cfg.loghYBlockSize;

            switch yGridMode
                case "fixed-bound"
                    cfg.AMI_Evaluator = @(x_in) AMI_noCSI_fast_logh_grid( ...
                        x_in, params.px, params, dt_local, span_local, ...
                        y_grid_local, block_local);

                case "candidate-adaptive"
                    scale_local = cfg.yGridScale;
                    cfg.AMI_Evaluator = @(x_in) AMI_noCSI_fast_logh_adaptive_grid( ...
                        x_in, params.px, params, dt_local, span_local, ...
                        scale_local, block_local);
            end
    end

    cfg.AMI_Validator = @(x_in) AMI_functions.AMI_noCSI_validate( ...
        x_in, params.px, params);
end
