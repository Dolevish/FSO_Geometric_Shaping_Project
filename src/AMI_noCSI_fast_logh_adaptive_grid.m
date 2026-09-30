function [mi_bits, diagnostics] = AMI_noCSI_fast_logh_adaptive_grid(x, px, params, dt, spanSigma, yGridScale, yBlockSize)
%AMI_NOCSI_FAST_LOGH_ADAPTIVE_GRID  Log-h AMI with candidate-adaptive y-grid.
%
%   mi_bits = AMI_noCSI_fast_logh_adaptive_grid(x,px,params,...)
%   [mi_bits,diagnostics] = ...
%
% Commit I.6 wrapper around AMI_noCSI_fast_logh_grid().  The fading
% quadrature is unchanged.  Only the outer y-grid support is rebuilt for the
% exact candidate constellation being scored:
%
%       x_support = yGridScale * max(abs(x)).
%
% The constellation itself is never clipped, projected, rescaled, or
% rejected by this function.  Therefore this mode does NOT introduce a peak
% power constraint.  Commit I.5 selected yGridScale=1.1 as the candidate for
% optimizer-level A/B validation.
%
% Commit I.10.1 performance path:
%   - channel-only y-grid constants are cached per MATLAB worker;
%   - diagnostics/timing are skipped completely for the scalar-output hot
%     path used by SA;
%   - the same support formulas, spacings and unique() construction are kept,
%     so the numerical y-grid is unchanged.
%
% Defaults:
%   dt          = 0.01
%   spanSigma   = 8
%   yGridScale  = 1.1
%   yBlockSize  = 512

    if nargin < 4 || isempty(dt), dt = 0.01; end
    if nargin < 5 || isempty(spanSigma), spanSigma = 8; end
    if nargin < 6 || isempty(yGridScale), yGridScale = 1.1; end
    if nargin < 7 || isempty(yBlockSize), yBlockSize = 512; end

    validateattributes(yGridScale,{'numeric'}, ...
        {'scalar','real','finite','>=',1},mfilename,'yGridScale');

    xBefore = real(x(:));
    xSupport = double(yGridScale) * max(abs(xBefore));
    if ~(isfinite(xSupport) && xSupport > 0)
        error('AMI_noCSI_fast_logh_adaptive_grid:BadSupport', ...
            'Candidate constellation must imply a finite positive y-grid support.');
    end

    if nargout > 1
        tGrid = tic;
        yGrid = cached_adaptive_y_grid(params,xSupport);
        gridBuildSeconds = toc(tGrid);
        [mi_bits,inner] = AMI_noCSI_fast_logh_grid( ...
            x,px,params,dt,spanSigma,yGrid,yBlockSize);

        % Defensive invariant is retained on the diagnostic path. The scalar
        % production path does not copy xBefore back through an additional
        % isequaln() check on every SA evaluation.
        if ~isequaln(real(x(:)),xBefore)
            error('AMI_noCSI_fast_logh_adaptive_grid:ConstellationModified', ...
                'Adaptive evaluator modified the candidate constellation.');
        end

        diagnostics = inner;
        diagnostics.yGridMode = 'candidate-adaptive';
        diagnostics.yGridScale = double(yGridScale);
        diagnostics.xMax = max(abs(xBefore));
        diagnostics.xSupport = xSupport;
        diagnostics.yGridN = numel(yGrid);
        diagnostics.yGridMin = yGrid(1);
        diagnostics.yGridMax = yGrid(end);
        diagnostics.gridBuildSeconds = gridBuildSeconds;
    else
        % SA hot path: no stopwatch, no diagnostics struct, no defensive
        % post-call copy comparison. The likelihood calculation is identical.
        yGrid = cached_adaptive_y_grid(params,xSupport);
        mi_bits = AMI_noCSI_fast_logh_grid( ...
            x,px,params,dt,spanSigma,yGrid,yBlockSize);
    end
end


function y_grid = cached_adaptive_y_grid(params,x_max_bound)
%CACHED_ADAPTIVE_Y_GRID  Same Commit-I.6 grid using cached channel constants.
% The old builder recomputed norminv() quantiles on every candidate. They
% depend only on the physical channel, so cache them once per worker/case.

    persistent lastKey lastPrep

    key = [double(params.sigma_n_sq),double(params.sig_t), ...
        double(params.mu_t),double(params.R)];

    if isempty(lastKey) || ~isequal(lastKey,key) || isempty(lastPrep)
        sig_n = sqrt(params.sigma_n_sq);
        if params.sig_t > 1e-6
            h_75   = exp(params.mu_t + norminv(0.75)   * params.sig_t);
            h_999  = exp(params.mu_t + norminv(0.999)  * params.sig_t);
            h_9999 = exp(params.mu_t + norminv(0.9999) * params.sig_t);
        else
            h_75 = 1; h_999 = 1; h_9999 = 1;
        end

        prep = struct();
        prep.sig_n = sig_n;
        prep.h_75 = h_75;
        prep.h_999 = h_999;
        prep.h_9999 = h_9999;
        prep.R = params.R;
        prep.y_min = -8*sig_n;
        prep.dy_fine = max(sig_n/4,1e-6);
        prep.dy_medium = max(sig_n,prep.dy_fine*4);
        prep.dy_coarse = max(sig_n*4,prep.dy_medium*4);

        lastKey = key;
        lastPrep = prep;
    else
        prep = lastPrep;
    end

    sig_n = prep.sig_n;
    y_focus = prep.R*prep.h_75*x_max_bound + 8*sig_n;
    y_tail  = prep.R*prep.h_999*x_max_bound + 5*sig_n;
    y_max   = prep.R*prep.h_9999*x_max_bound + 10*sig_n;

    y_focus = max(y_focus,prep.y_min + 20*sig_n);
    y_tail  = max(y_tail,y_focus + 10*sig_n);
    y_max   = max(y_max,y_tail + 10*sig_n);

    g_dense  = prep.y_min : prep.dy_fine   : y_focus;
    g_medium = (y_focus + prep.dy_medium) : prep.dy_medium : y_tail;
    g_coarse = (y_tail  + prep.dy_coarse) : prep.dy_coarse : y_max;

    % Keep the historical unique() operation intentionally. This commit is a
    % performance refactor, not a numerical-grid change.
    y_grid = unique([g_dense,g_medium,g_coarse]);
end
