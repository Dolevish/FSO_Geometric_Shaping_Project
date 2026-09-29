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

    tGrid = tic;
    yGrid = AMI_functions.build_noCSI_y_grid(params,xSupport);
    gridBuildSeconds = toc(tGrid);

    [mi_bits,inner] = AMI_noCSI_fast_logh_grid( ...
        x,px,params,dt,spanSigma,yGrid,yBlockSize);

    % Defensive invariant: evaluator must be side-effect free.
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
end
