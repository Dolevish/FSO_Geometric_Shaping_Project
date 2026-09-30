function [mi_bits, diagnostics] = AMI_noCSI_fast_logh_adaptive_grid(x, px, params, dt, spanSigma, yGridScale, yBlockSize)
%AMI_NOCSI_FAST_LOGH_ADAPTIVE_GRID  Log-h AMI with candidate-adaptive y-grid.
%
% Commit I.10.2 runtime path keeps the exact I.6 support formulas while:
%   - using deterministic M-dependent y-block tuning by default;
%   - bypassing unique() because the dense/medium/coarse grid segments are
%     constructed strictly separated and already increasing;
%   - delegating likelihood work to the I.10.2 kernel.
%
% Set params.loghRuntimeKernel='i10p1' to reproduce the I.10.1 path.
% Set params.loghRuntimeBlockPolicy='fixed' to force the supplied yBlockSize.
% The default block policy for the I.10.2 kernel is 'auto'.
%
% None of these implementation choices changes the candidate support,
% integration spacings, or constellation being scored.

    if nargin < 4 || isempty(dt), dt = 0.01; end
    if nargin < 5 || isempty(spanSigma), spanSigma = 8; end
    if nargin < 6 || isempty(yGridScale), yGridScale = 1.1; end
    if nargin < 7 || isempty(yBlockSize), yBlockSize = 512; end

    validateattributes(yGridScale,{'numeric'}, ...
        {'scalar','real','finite','>=',1},mfilename,'yGridScale');
    validateattributes(yBlockSize,{'numeric'}, ...
        {'scalar','real','finite','integer','positive'},mfilename,'yBlockSize');

    xBefore = real(x(:));
    xSupport = double(yGridScale) * max(abs(xBefore));
    if ~(isfinite(xSupport) && xSupport > 0)
        error('AMI_noCSI_fast_logh_adaptive_grid:BadSupport', ...
            'Candidate constellation must imply a finite positive y-grid support.');
    end

    kernelMode=runtime_kernel_mode(params);
    [effectiveBlock,blockPolicy,tuningMeta]=resolve_block_size(params,numel(xBefore),yBlockSize,kernelMode);
    fastGrid = kernelMode=="i10p2";

    if nargout > 1
        tGrid = tic;
        yGrid = cached_adaptive_y_grid(params,xSupport,fastGrid);
        gridBuildSeconds = toc(tGrid);
        [mi_bits,inner] = AMI_noCSI_fast_logh_grid( ...
            x,px,params,dt,spanSigma,yGrid,effectiveBlock);

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
        diagnostics.yBlockSizeRequested = double(yBlockSize);
        diagnostics.yBlockSizeEffective = double(effectiveBlock);
        diagnostics.runtimeBlockPolicy = char(blockPolicy);
        diagnostics.runtimeTuningRuleVersion = tuningMeta.ruleVersion;
        diagnostics.fastGridConcat = logical(fastGrid);
    else
        yGrid = cached_adaptive_y_grid(params,xSupport,fastGrid);
        mi_bits = AMI_noCSI_fast_logh_grid( ...
            x,px,params,dt,spanSigma,yGrid,effectiveBlock);
    end
end


function mode=runtime_kernel_mode(params)
    mode="i10p2";
    if isfield(params,'loghRuntimeKernel') && ~isempty(params.loghRuntimeKernel)
        mode=lower(string(params.loghRuntimeKernel));
    end
    if ~ismember(mode,["i10p1","i10p2"])
        error('AMI_noCSI_fast_logh_adaptive_grid:BadRuntimeKernel', ...
            'loghRuntimeKernel must be ''i10p1'' or ''i10p2''.');
    end
end


function [effectiveBlock,policy,meta]=resolve_block_size(params,M,requestedBlock,kernelMode)
    [autoBlock,~,meta]=select_logh_runtime_tuning(M);
    policy="auto";

    if kernelMode=="i10p1"
        policy="fixed";
    elseif isfield(params,'loghRuntimeBlockPolicy') && ~isempty(params.loghRuntimeBlockPolicy)
        policy=lower(string(params.loghRuntimeBlockPolicy));
    end

    if ~ismember(policy,["auto","fixed"])
        error('AMI_noCSI_fast_logh_adaptive_grid:BadRuntimeBlockPolicy', ...
            'loghRuntimeBlockPolicy must be ''auto'' or ''fixed''.');
    end

    if policy=="auto"
        effectiveBlock=autoBlock;
    else
        effectiveBlock=double(requestedBlock);
    end
end


function y_grid = cached_adaptive_y_grid(params,x_max_bound,fastConcat)
%CACHED_ADAPTIVE_Y_GRID  Same physical y-grid with cached channel constants.
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

    g_dense  = prep.y_min : prep.dy_fine : y_focus;
    g_medium = (y_focus + prep.dy_medium) : prep.dy_medium : y_tail;
    g_coarse = (y_tail + prep.dy_coarse) : prep.dy_coarse : y_max;

    if fastConcat
        % Each later segment starts at least one positive step beyond the
        % previous support boundary, so the concatenation is already strictly
        % increasing and contains no duplicates. Therefore unique() was pure
        % overhead in the SA hot path.
        y_grid=[g_dense,g_medium,g_coarse];
    else
        % Exact historical I.10.1 path retained for regression benchmarks.
        y_grid=unique([g_dense,g_medium,g_coarse]);
    end
end
