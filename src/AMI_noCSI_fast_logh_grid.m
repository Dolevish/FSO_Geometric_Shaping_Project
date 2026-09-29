function [mi_bits, diagnostics] = AMI_noCSI_fast_logh_grid(x, px, params, dt, spanSigma, y_grid, yBlockSize)
%AMI_NOCSI_FAST_LOGH_GRID  Fast no-CSI AMI via uniform quadrature in t=ln(h).
%
%   mi_bits = AMI_noCSI_fast_logh_grid(x,px,params)
%   mi_bits = AMI_noCSI_fast_logh_grid(x,px,params,dt,spanSigma,y_grid,yBlockSize)
%   [mi_bits,diagnostics] = ...
%
% Channel model:
%       y = h R x + n,
%       t = ln(h) ~ N(mu_t,sig_t^2),
%       n ~ N(0,sigma_n^2).
%
% The fading mixture is evaluated directly on a uniform grid in t=ln(h):
%
%   p(y|x) = integral p(y|x,t) f_T(t) dt,
%
% over [mu_t-spanSigma*sig_t, mu_t+spanSigma*sig_t].  Trapezoidal
% node weights are used in t.  The outer AMI integral uses the supplied
% non-uniform y-grid and trapz, matching the project's fast evaluator.
%
% Defaults selected by Commit I.2 diagnostics:
%   dt         = 0.01
%   spanSigma  = 8
%   yBlockSize = 512
%
% The y dimension is processed in blocks to bound temporary memory usage.
% This function is side-effect free and scores the exact constellation x.

    if nargin < 4 || isempty(dt), dt = 0.01; end
    if nargin < 5 || isempty(spanSigma), spanSigma = 8; end
    if nargin < 6 || isempty(y_grid)
        y_grid = AMI_functions.build_noCSI_y_grid(params, max(abs(x))*1.5);
    end
    if nargin < 7 || isempty(yBlockSize), yBlockSize = 512; end

    validateattributes(dt,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(spanSigma,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(yBlockSize,{'numeric'},{'scalar','integer','positive'});

    x = real(x(:).');
    px = real(px(:).');
    y_grid = real(y_grid(:).');

    M = numel(x);
    Ny = numel(y_grid);

    if isempty(x) || any(~isfinite(x))
        error('AMI_noCSI_fast_logh_grid:BadConstellation', ...
            'x must be non-empty and contain finite values only.');
    end
    if numel(px) ~= M || any(~isfinite(px)) || any(px < 0) || sum(px) <= 0
        error('AMI_noCSI_fast_logh_grid:BadProbabilities', ...
            'px must contain one finite nonnegative probability per constellation point.');
    end
    px = px / sum(px);
    if numel(y_grid) < 2 || any(~isfinite(y_grid)) || any(diff(y_grid) <= 0)
        error('AMI_noCSI_fast_logh_grid:BadYGrid', ...
            'y_grid must be a finite strictly increasing vector with at least two points.');
    end

    % Degenerate no-turbulence case: h=1, so the fading integral disappears.
    if params.sig_t < 1e-12
        means = params.R * x(:);
        d2 = bsxfun(@minus,y_grid,means).^2;
        logPyx = -0.5*log(2*pi*params.sigma_n_sq) - d2/(2*params.sigma_n_sq);
        mi_bits = mi_from_log_pyx(logPyx,px,y_grid);
        diagnostics = struct('nT',1,'capturedMass',1,'maxStep',0, ...
            'tMin',params.mu_t,'tMax',params.mu_t,'dtRequested',double(dt), ...
            'spanSigma',double(spanSigma),'yBlockSize',double(yBlockSize));
        return;
    end

    tMin = params.mu_t - spanSigma*params.sig_t;
    tMax = params.mu_t + spanSigma*params.sig_t;

    tGrid = tMin:dt:tMax;
    if isempty(tGrid) || tGrid(end) < tMax - 100*eps(max(1,abs(tMax)))
        tGrid = [tGrid tMax]; %#ok<AGROW>
    else
        tGrid(end) = tMax;
    end
    tGrid = unique(tGrid,'stable');

    qWeights = trapezoid_node_weights(tGrid);
    fT = exp(-0.5*((tGrid-params.mu_t)/params.sig_t).^2) / ...
        (sqrt(2*pi)*params.sig_t);
    q = qWeights .* fT;
    capturedMass = sum(q);

    keep = q > 0 & isfinite(q);
    tGrid = tGrid(keep);
    q = q(keep);

    if isempty(tGrid) || ~(capturedMass > 0)
        error('AMI_noCSI_fast_logh_grid:EmptyQuadrature', ...
            'Uniform log-h quadrature produced no positive finite weights.');
    end

    logQ = log(q(:));
    hVals = exp(tGrid(:));

    inv2s2 = 1/(2*params.sigma_n_sq);
    logNorm = -0.5*log(2*pi*params.sigma_n_sq);
    logPyx = zeros(M,Ny);

    for j = 1:M
        means = params.R*hVals*x(j);
        for b0 = 1:yBlockSize:Ny
            b1 = min(Ny,b0+yBlockSize-1);
            yb = y_grid(b0:b1);
            d2 = bsxfun(@minus,yb,means).^2;
            logTerms = bsxfun(@plus,logQ+logNorm,-inv2s2*d2);
            logPyx(j,b0:b1) = local_logsumexp(logTerms,1);
        end
    end

    mi_bits = mi_from_log_pyx(logPyx,px,y_grid);

    steps = diff(tGrid);
    if isempty(steps), maxStep = 0; else, maxStep = max(steps); end
    diagnostics = struct('nT',numel(tGrid),'capturedMass',capturedMass, ...
        'maxStep',maxStep,'tMin',tMin,'tMax',tMax, ...
        'dtRequested',double(dt),'spanSigma',double(spanSigma), ...
        'yBlockSize',double(yBlockSize));
end


function miBits = mi_from_log_pyx(logPyx,px,yGrid)
    logPxPyx = bsxfun(@plus,log(px(:)),logPyx);
    logPy = local_logsumexp(logPxPyx,1);
    M = size(logPyx,1);
    MI = 0;
    for j = 1:M
        Pyxj = exp(logPyx(j,:));
        lr = (logPyx(j,:)-logPy)/log(2);
        integrand = Pyxj.*lr;
        integrand(~isfinite(integrand)) = 0;
        MI = MI + px(j)*trapz(yGrid,integrand);
    end
    miBits = max(0,MI);
end


function w = trapezoid_node_weights(x)
    x = x(:).';
    n = numel(x);
    if n == 1
        w = 1;
        return;
    end
    dx = diff(x);
    w = zeros(1,n);
    w(1) = dx(1)/2;
    w(end) = dx(end)/2;
    if n > 2
        w(2:end-1) = (dx(1:end-1)+dx(2:end))/2;
    end
end


function s = local_logsumexp(A,dim)
    if nargin < 2, dim = 1; end
    m = max(A,[],dim);
    s = m + log(sum(exp(bsxfun(@minus,A,m)),dim));
end
