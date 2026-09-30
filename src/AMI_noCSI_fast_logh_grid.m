function [mi_bits, diagnostics] = AMI_noCSI_fast_logh_grid(x, px, params, dt, spanSigma, y_grid, yBlockSize)
%AMI_NOCSI_FAST_LOGH_GRID  Fast no-CSI AMI via uniform quadrature in t=ln(h).
%
% Channel model:
%       y = h R x + n,
%       t = ln(h) ~ N(mu_t,sig_t^2),
%       n ~ N(0,sigma_n^2).
%
% Commit I.10.1 caches all channel-dependent log-h quadrature constants.
%
% Commit I.10.2 optimizes only the likelihood kernel:
%   1) exact x=0 fast path (including the finite quadrature captured mass),
%   2) small memory-aware symbol batches,
%   3) lower-allocation work arrays in the hot path.
%
% The historical I.10.1 kernel remains available for regression/benchmark by
% setting params.loghRuntimeKernel='i10p1'.  The default is 'i10p2'.
% Optional implementation-only overrides:
%   params.loghSymbolBatchSize   positive integer
%   params.loghZeroFastPath      logical scalar
%
% None of these settings changes the t-grid, y-grid, probability model,
% objective definition, or constellation itself.

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

    [kernelMode,symbolBatchSize,zeroFastPath] = runtime_kernel_settings(params,M);

    % Degenerate no-turbulence case: h=1, so the fading integral disappears.
    if params.sig_t < 1e-12
        means = params.R * x(:);
        d2 = bsxfun(@minus,y_grid,means).^2;
        logPyx = -0.5*log(2*pi*params.sigma_n_sq) - d2/(2*params.sigma_n_sq);
        mi_bits = mi_from_log_pyx(logPyx,px,y_grid);
        if nargout > 1
            diagnostics = struct('nT',1,'capturedMass',1,'maxStep',0, ...
                'tMin',params.mu_t,'tMax',params.mu_t,'dtRequested',double(dt), ...
                'spanSigma',double(spanSigma),'yBlockSize',double(yBlockSize), ...
                'quadratureCacheHit',false,'kernelMode',char(kernelMode), ...
                'symbolBatchSize',double(symbolBatchSize),'zeroFastPath',logical(zeroFastPath), ...
                'nZeroFastSymbols',0);
        end
        return;
    end

    [prep,cacheHit] = cached_logh_quadrature(params,dt,spanSigma);

    if kernelMode == "i10p1"
        [logPyx,nZeroFast] = likelihood_i10p1(x,prep,y_grid,yBlockSize);
    else
        [logPyx,nZeroFast] = likelihood_i10p2(x,prep,y_grid,yBlockSize, ...
            symbolBatchSize,zeroFastPath);
    end

    mi_bits = mi_from_log_pyx(logPyx,px,y_grid);

    if nargout > 1
        diagnostics = struct('nT',prep.nT,'capturedMass',prep.capturedMass, ...
            'maxStep',prep.maxStep,'tMin',prep.tMin,'tMax',prep.tMax, ...
            'dtRequested',double(dt),'spanSigma',double(spanSigma), ...
            'yBlockSize',double(yBlockSize),'quadratureCacheHit',cacheHit, ...
            'kernelMode',char(kernelMode),'symbolBatchSize',double(symbolBatchSize), ...
            'zeroFastPath',logical(zeroFastPath),'nZeroFastSymbols',double(nZeroFast));
    end
end


function [logPyx,nZeroFast] = likelihood_i10p1(x,prep,yGrid,yBlockSize)
% Exact Commit-I.10.1 per-symbol/per-block likelihood loop.
    M=numel(x); Ny=numel(yGrid); logPyx=zeros(M,Ny); nZeroFast=0;
    for j=1:M
        means=prep.hScale*x(j);
        for b0=1:yBlockSize:Ny
            b1=min(Ny,b0+yBlockSize-1); yb=yGrid(b0:b1);
            d2=bsxfun(@minus,yb,means).^2;
            logTerms=bsxfun(@plus,prep.logQNorm,-prep.inv2s2*d2);
            logPyx(j,b0:b1)=local_logsumexp(logTerms,1);
        end
    end
end


function [logPyx,nZeroFast] = likelihood_i10p2(x,prep,yGrid,yBlockSize,batchSize,zeroFastPath)
% Memory-aware hot kernel. Arithmetic model and quadrature are unchanged.
    M=numel(x); Ny=numel(yGrid); logPyx=zeros(M,Ny);

    zeroIdx=[];
    if zeroFastPath
        zeroIdx=find(x==0);
        if ~isempty(zeroIdx)
            % Historical kernel for x=0 is LSE(logQNorm - inv2s2*y^2).
            % Since the y term is common to all t nodes, factor it out while
            % preserving the finite quadrature captured mass through the same
            % cached log-sum-exp normalization.
            zeroLog = prep.logIntegratedNorm - prep.inv2s2*(yGrid.^2);
            logPyx(zeroIdx,:) = repmat(zeroLog,numel(zeroIdx),1);
        end
    end
    nZeroFast=numel(zeroIdx);

    active=true(1,M); active(zeroIdx)=false; activeIdx=find(active);
    if isempty(activeIdx), return; end

    batchSize=max(1,min(round(double(batchSize)),numel(activeIdx)));

    % y-block outer loop keeps each temporary bounded.  A batch of symbols is
    % flattened into adjacent column groups [y_1..y_B] for each x_j, allowing
    % MATLAB to amortize loop/logsumexp overhead without a large 3-D array.
    for b0=1:yBlockSize:Ny
        b1=min(Ny,b0+yBlockSize-1);
        yb=yGrid(b0:b1); nb=numel(yb);

        for a0=1:batchSize:numel(activeIdx)
            idx=activeIdx(a0:min(numel(activeIdx),a0+batchSize-1));
            nS=numel(idx);

            if nS==1
                means=prep.hScale*x(idx);
                work=bsxfun(@minus,yb,means);
                work=work.^2;
                work=(-prep.inv2s2).*work;
                work=bsxfun(@plus,prep.logQNorm,work);
                logPyx(idx,b0:b1)=local_logsumexp(work,1);
            else
                yCols=repmat(yb,1,nS);
                xCols=repelem(x(idx),nb);
                work=bsxfun(@times,prep.hScale,xCols);
                work=bsxfun(@minus,yCols,work);
                work=work.^2;
                work=(-prep.inv2s2).*work;
                work=bsxfun(@plus,prep.logQNorm,work);
                L=local_logsumexp(work,1);
                logPyx(idx,b0:b1)=reshape(L,nb,nS).';
            end
        end
    end
end


function [mode,batchSize,zeroFastPath] = runtime_kernel_settings(params,M)
    mode="i10p2";
    if isfield(params,'loghRuntimeKernel') && ~isempty(params.loghRuntimeKernel)
        mode=lower(string(params.loghRuntimeKernel));
    end
    if ~ismember(mode,["i10p1","i10p2"])
        error('AMI_noCSI_fast_logh_grid:BadRuntimeKernel', ...
            'loghRuntimeKernel must be ''i10p1'' or ''i10p2''.');
    end

    [~,autoBatch]=select_logh_runtime_tuning(M);
    batchSize=autoBatch;
    zeroFastPath=true;

    if mode=="i10p1"
        batchSize=1; zeroFastPath=false;
        return;
    end
    if isfield(params,'loghSymbolBatchSize') && ~isempty(params.loghSymbolBatchSize)
        validateattributes(params.loghSymbolBatchSize,{'numeric'}, ...
            {'scalar','real','finite','integer','positive'});
        batchSize=double(params.loghSymbolBatchSize);
    end
    batchSize=max(1,min(batchSize,M));

    if isfield(params,'loghZeroFastPath') && ~isempty(params.loghZeroFastPath)
        if ~(islogical(params.loghZeroFastPath) && isscalar(params.loghZeroFastPath))
            error('AMI_noCSI_fast_logh_grid:BadZeroFastPath', ...
                'loghZeroFastPath must be a logical scalar.');
        end
        zeroFastPath=params.loghZeroFastPath;
    end
end


function [prep,cacheHit] = cached_logh_quadrature(params,dt,spanSigma)
%CACHED_LOGH_QUADRATURE  One-entry per-worker cache for a physical channel case.
    persistent lastKey lastPrep

    key = [double(params.mu_t),double(params.sig_t),double(params.sigma_n_sq), ...
        double(params.R),double(dt),double(spanSigma)];

    cacheHit = ~isempty(lastKey) && isequal(lastKey,key) && ~isempty(lastPrep);
    if cacheHit
        prep = lastPrep;
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

    logNorm = -0.5*log(2*pi*params.sigma_n_sq);
    steps = diff(tGrid);
    if isempty(steps), maxStep = 0; else, maxStep = max(steps); end

    prep = struct();
    prep.logQNorm = log(q(:)) + logNorm;
    prep.logIntegratedNorm = local_logsumexp(prep.logQNorm,1);
    prep.hScale = params.R*exp(tGrid(:));
    prep.inv2s2 = 1/(2*params.sigma_n_sq);
    prep.capturedMass = capturedMass;
    prep.nT = numel(tGrid);
    prep.maxStep = maxStep;
    prep.tMin = tMin;
    prep.tMax = tMax;

    lastKey = key;
    lastPrep = prep;
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
    if n == 1, w = 1; return; end
    dx = diff(x);
    w = zeros(1,n);
    w(1) = dx(1)/2;
    w(end) = dx(end)/2;
    if n > 2, w(2:end-1) = (dx(1:end-1)+dx(2:end))/2; end
end


function s = local_logsumexp(A,dim)
    if nargin < 2, dim = 1; end
    m = max(A,[],dim);
    s = m + log(sum(exp(bsxfun(@minus,A,m)),dim));
end
