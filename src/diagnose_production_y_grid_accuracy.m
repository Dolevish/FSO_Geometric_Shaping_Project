function bundle = diagnose_production_y_grid_accuracy(varargin)
%DIAGNOSE_PRODUCTION_Y_GRID_ACCURACY  Commit-I.7 full-grid no-SA accuracy sweep.
%
% Complements the optimizer-level Commit-I.6 A/B test by checking the
% candidate-adaptive production y-grid across the entire manuscript grid,
% without rerunning SA.  Each physical point is evaluated on deterministic
% feasible constellations that span moderate and strongly clustered shapes.
%
% Production defaults:
%   M                 = [4 8 16 32]
%   SNR               = 5:5:30 dB
%   sigma_X^2          = [0 0.1 0.2 0.3]
%   P_avg             = 1
%   d_min             = 0.01
%   log-h dt          = 0.01
%   log-h span        = mu_t +/- 8 sigma_t
%   adaptive scale    = 1.1
%   error tolerance   = 1e-3 bit/symbol
%
% Deterministic geometry set:
%   uniform        : ordinary positive Uniform-PAM with x_1=0, mean=P_avg.
%   upper-quarter  : minimum-gap skeleton plus all surplus power shared by
%                    the largest ceil(M/4) levels.
%   upper-eighth   : same construction using ceil(M/8) largest levels.
%
% The clustered geometries preserve x_1=0, sorting, d_min, and mean power
% exactly. Duplicate geometries for small M are removed per physical case.
% No clipping, projection, or peak constraint is introduced.
%
% The independent validator is the existing adaptive-integral validator.
% Both fixed-bound and candidate-adaptive log-h fast evaluators are compared
% against exactly the same validator result for every geometry.
%
% This diagnostic is intentionally NO-SA: it is a pre-production numerical
% coverage check, not an optimization-quality experiment.

    p = inputParser;
    p.FunctionName = mfilename;
    addParameter(p,'MVec',[4 8 16 32],@positive_integer_vector);
    addParameter(p,'P_avg',1,@positive_scalar);
    addParameter(p,'SNRVec',5:5:30,@finite_vector);
    addParameter(p,'TurbulenceVec',[0 0.1 0.2 0.3],@nonnegative_vector);
    addParameter(p,'minGap',0.01,@nonnegative_scalar);
    addParameter(p,'pinZero',true,@(v) islogical(v)&&isscalar(v));
    addParameter(p,'loghDt',0.01,@positive_scalar);
    addParameter(p,'loghSpanSigma',8,@positive_scalar);
    addParameter(p,'loghYBlockSize',512,@positive_integer);
    addParameter(p,'AdaptiveScale',1.1,@(v) positive_scalar(v)&&v>=1);
    addParameter(p,'AbsErrorTolerance',1e-3,@positive_scalar);
    addParameter(p,'GeometryModes',{'uniform','upper-quarter','upper-eighth'},@geometry_modes_valid);
    addParameter(p,'UseParallelCases',true,@(v) islogical(v)&&isscalar(v));
    addParameter(p,'WriteCSV',true,@(v) islogical(v)&&isscalar(v));
    addParameter(p,'RunLabel','',@(v) ischar(v)||(isstring(v)&&isscalar(v)));
    addParameter(p,'ResultsRoot',fso_result_utils.default_results_root(),@(v) ischar(v)||isstring(v));
    parse(p,varargin{:});
    o = p.Results;

    if ~o.pinZero
        error('diagnose_production_y_grid_accuracy:PinZeroRequired', ...
            'Commit-I.7 deterministic stress geometries are defined for pinZero=true.');
    end

    MVec = unique(double(o.MVec(:).'),'stable');
    snrVec = unique(double(o.SNRVec(:).'),'stable');
    sigVec = unique(double(o.TurbulenceVec(:).'),'stable');
    modes = cellstr(string(o.GeometryModes(:).'));
    Pavg = double(o.P_avg);
    dmin = double(o.minGap);
    validate_spacing(MVec,Pavg,dmin);

    [iM,iSig,iSNR] = ndgrid(1:numel(MVec),1:numel(sigVec),1:numel(snrVec));
    tM = iM(:); tSig = iSig(:); tSNR = iSNR(:);
    nPhysical = numel(tM);

    resultsRoot = char(o.ResultsRoot);
    parentDir = fullfile(resultsRoot,'production_y_grid_accuracy');
    label = sanitize_label(o.RunLabel);
    token = sprintf('LOGHdt%s_T%s_scale%s_tol%s', ...
        num_token(o.loghDt,4),num_token(o.loghSpanSigma,2), ...
        num_token(o.AdaptiveScale,2),num_token(o.AbsErrorTolerance,5));
    utcToken = char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmssSSS''Z'''));
    [outDir,runId] = make_unique_dir(parentDir,sprintf('%s_%s_%s',label,token,utcToken));

    runMeta = fso_result_utils.run_metadata(mfilename);
    runMeta.runId = runId;
    runMeta.runLabel = char(string(o.RunLabel));
    runMeta.outputDirectory = outDir;

    fprintf('\n============================================================\n');
    fprintf('Commit I.7 - production candidate-adaptive y-grid accuracy sweep\n');
    fprintf('Run ID: %s\n',runId);
    fprintf('M=%s | SNR=%s dB | sigma_X^2=%s | Pavg=%.4g | d_min=%.4g\n', ...
        mat2str(MVec),mat2str(snrVec),mat2str(sigVec),Pavg,dmin);
    fprintf('log-h dt=%.5g | span=+/-%.3g sigma_t | yBlock=%d\n', ...
        o.loghDt,o.loghSpanSigma,o.loghYBlockSize);
    fprintf('Adaptive scale=%.4g | tolerance=%.3e bit/symbol\n', ...
        o.AdaptiveScale,o.AbsErrorTolerance);
    fprintf('Geometry modes=%s\n',strjoin(modes,', '));
    fprintf('Physical points=%d | no SA is run\n',nPhysical);
    fprintf('Results: %s\n',outDir);
    fprintf('============================================================\n\n');

    caseRows = cell(nPhysical,1);
    doPar = logical(o.UseParallelCases);
    if doPar
        try
            doPar = license('test','Distrib_Computing_Toolbox') ~= 0;
        catch
            doPar = false;
        end
    end

    if doPar
        pool = gcp('nocreate');
        if isempty(pool), pool = parpool('Processes'); end
        fprintf('Outer physical-point parallelism: %d process workers\n',pool.NumWorkers);
        Q = parallel.pool.DataQueue;
        doneCount = 0;
        afterEach(Q,@(~) tick());
        parfor k = 1:nPhysical
            caseRows{k} = evaluate_point(MVec(tM(k)),snrVec(tSNR(k)),sigVec(tSig(k)),Pavg,dmin,modes,o);
            send(Q,1);
        end
    else
        fprintf('Outer physical-point parallelism: disabled\n');
        for k = 1:nPhysical
            caseRows{k} = evaluate_point(MVec(tM(k)),snrVec(tSNR(k)),sigVec(tSig(k)),Pavg,dmin,modes,o);
            fprintf('  completed %d/%d\n',k,nPhysical);
        end
    end

    T = vertcat(caseRows{:});
    summary = summarize_table(T,double(o.AbsErrorTolerance));

    bundle = struct();
    bundle.meta = runMeta;
    bundle.experiment = 'production_y_grid_accuracy';
    bundle.run = struct('id',runId,'label',char(string(o.RunLabel)), ...
        'outputDirectory',outDir);
    bundle.axes = struct('M',MVec,'P_avg',Pavg,'SNR_dB',snrVec, ...
        'sigma_X_sq',sigVec,'geometryModes',{modes});
    bundle.method = struct('fastFadingMethod','logh','loghDt',double(o.loghDt), ...
        'loghSpanSigma',double(o.loghSpanSigma),'loghYBlockSize',double(o.loghYBlockSize), ...
        'fixedYGridMode','fixed-bound','adaptiveYGridMode','candidate-adaptive', ...
        'adaptiveScale',double(o.AdaptiveScale),'absErrorTolerance',double(o.AbsErrorTolerance));
    bundle.results = T;
    bundle.summary = summary;

    matPath = fullfile(outDir,'productionYGridAccuracy.mat');
    save(matPath,'bundle','-v7.3');
    if o.WriteCSV
        writetable(T,fullfile(outDir,'productionYGridAccuracy.csv'));
    end

    fprintf('\nSaved I.7 accuracy bundle: %s\n',matPath);
    print_summary(summary,T);

    function tick()
        doneCount = doneCount + 1;
        fprintf('  completed %d/%d\n',doneCount,nPhysical);
    end
end


function T = evaluate_point(M,SNRdB,sigmaX2,Pavg,dmin,modes,o)
    xMaxBound = fso_result_utils.feasible_xmax_bound(M,Pavg,dmin,true);
    common = {'FastFadingMethod','logh','loghDt',o.loghDt, ...
        'loghSpanSigma',o.loghSpanSigma,'loghYBlockSize',o.loghYBlockSize, ...
        'xMaxBound',xMaxBound,'saMaxIter',1,'saNStarts',1,'saUseParallel',false, ...
        'minGap',dmin,'pinZero',true,'seedInit',1,'logEvery',0,'historyEvery',1};

    cfgFixed = build_fso_config(M,Pavg,SNRdB,sigmaX2,common{:}, ...
        'YGridMode','fixed-bound','YGridScale',o.AdaptiveScale);
    cfgAdapt = build_fso_config(M,Pavg,SNRdB,sigmaX2,common{:}, ...
        'YGridMode','candidate-adaptive','YGridScale',o.AdaptiveScale);

    [X,names] = deterministic_geometries(M,Pavg,dmin,modes);
    nG = numel(names);

    Method = strings(nG,1); %#ok<NASGU>
    Geometry = strings(nG,1);
    Mcol = repmat(M,nG,1);
    SigmaX2 = repmat(sigmaX2,nG,1);
    SNRdBcol = repmat(SNRdB,nG,1);
    XMax = nan(nG,1);
    XMaxToBound = nan(nG,1);
    ValidatorAMI = nan(nG,1);
    FixedAMI = nan(nG,1);
    AdaptiveAMI = nan(nG,1);
    FixedError = nan(nG,1);
    AdaptiveError = nan(nG,1);
    AdaptiveMinusFixed = nan(nG,1);
    ValidatorSeconds = nan(nG,1);
    FixedSeconds = nan(nG,1);
    AdaptiveSeconds = nan(nG,1);
    FixedYGridN = repmat(numel(cfgFixed.y_grid),nG,1);
    AdaptiveYGridN = nan(nG,1);
    AdaptiveGridRatio = nan(nG,1);
    AdaptiveSupport = nan(nG,1);

    for g = 1:nG
        x = X{g};
        Geometry(g) = string(names{g});
        AMI_functions.assert_constellation_feasible(x,cfgFixed,1e-10);
        AMI_functions.assert_constellation_feasible(x,cfgAdapt,1e-10);
        XMax(g) = max(abs(x));
        XMaxToBound(g) = XMax(g)/xMaxBound;

        tv = tic;
        ValidatorAMI(g) = cfgFixed.AMI_Validator(x(:));
        ValidatorSeconds(g) = toc(tv);

        tf = tic;
        [FixedAMI(g),~] = AMI_noCSI_fast_logh_grid( ...
            x,cfgFixed.px,cfgFixed,o.loghDt,o.loghSpanSigma,cfgFixed.y_grid,o.loghYBlockSize);
        FixedSeconds(g) = toc(tf);

        ta = tic;
        [AdaptiveAMI(g),dA] = AMI_noCSI_fast_logh_adaptive_grid( ...
            x,cfgAdapt.px,cfgAdapt,o.loghDt,o.loghSpanSigma,o.AdaptiveScale,o.loghYBlockSize);
        AdaptiveSeconds(g) = toc(ta);

        FixedError(g) = FixedAMI(g)-ValidatorAMI(g);
        AdaptiveError(g) = AdaptiveAMI(g)-ValidatorAMI(g);
        AdaptiveMinusFixed(g) = AdaptiveAMI(g)-FixedAMI(g);
        AdaptiveYGridN(g) = dA.yGridN;
        AdaptiveGridRatio(g) = dA.yGridN/numel(cfgFixed.y_grid);
        AdaptiveSupport(g) = dA.xSupport;
    end

    T = table(Mcol,SigmaX2,SNRdBcol,Geometry,XMax,XMaxToBound, ...
        ValidatorAMI,FixedAMI,AdaptiveAMI,FixedError,AdaptiveError,AdaptiveMinusFixed, ...
        ValidatorSeconds,FixedSeconds,AdaptiveSeconds,FixedYGridN,AdaptiveYGridN, ...
        AdaptiveGridRatio,AdaptiveSupport, ...
        'VariableNames',{'M','SigmaX2','SNRdB','Geometry','XMax','XMaxToBound', ...
        'ValidatorAMI','FixedAMI','AdaptiveAMI','FixedError','AdaptiveError', ...
        'AdaptiveMinusFixed','ValidatorSeconds','FixedSeconds','AdaptiveSeconds', ...
        'FixedYGridN','AdaptiveYGridN','AdaptiveGridRatio','AdaptiveSupport'});
end


function [X,names] = deterministic_geometries(M,Pavg,dmin,modes)
    X = {};
    names = {};
    for j = 1:numel(modes)
        mode = lower(string(modes{j}));
        switch mode
            case "uniform"
                x = define_constellation(M,Pavg,true,"mean");
            case "upper-quarter"
                x = clustered_geometry(M,Pavg,dmin,max(1,ceil(M/4)));
            case "upper-eighth"
                x = clustered_geometry(M,Pavg,dmin,max(1,ceil(M/8)));
            case "upper-one"
                x = clustered_geometry(M,Pavg,dmin,1);
            otherwise
                error('diagnose_production_y_grid_accuracy:BadGeometryMode', ...
                    'Unknown geometry mode: %s',mode);
        end
        x = real(x(:));

        duplicate = false;
        for k = 1:numel(X)
            if numel(X{k})==numel(x) && max(abs(X{k}-x)) <= 1e-12*max(1,max(abs(x)))
                duplicate = true;
                break;
            end
        end
        if ~duplicate
            X{end+1,1} = x; %#ok<AGROW>
            names{end+1,1} = char(mode); %#ok<AGROW>
        end
    end
end


function x = clustered_geometry(M,Pavg,dmin,K)
    K = min(max(1,round(K)),M-1);
    b = (0:M-1).' * dmin;
    surplusMean = Pavg-mean(b);
    if surplusMean < -100*eps(max(1,Pavg))
        error('diagnose_production_y_grid_accuracy:InfeasibleGap', ...
            'Minimum-gap skeleton exceeds the requested average power.');
    end
    surplusMean = max(0,surplusMean);
    q = zeros(M,1);
    q(end-K+1:end) = M*surplusMean/K;
    x = b+q;
    x(1) = 0;
    % Correct only floating-point roundoff while preserving the intended form.
    err = Pavg-mean(x);
    if abs(err)>0
        x(end-K+1:end) = x(end-K+1:end) + M*err/K;
    end
end


function S = summarize_table(T,tol)
    absF = abs(T.FixedError);
    absA = abs(T.AdaptiveError);
    [S.maxAbsFixedError,idxF] = max(absF,[],'omitnan');
    [S.maxAbsAdaptiveError,idxA] = max(absA,[],'omitnan');
    S.meanAbsFixedError = mean(absF,'omitnan');
    S.meanAbsAdaptiveError = mean(absA,'omitnan');
    S.maxAbsAdaptiveMinusFixed = max(abs(T.AdaptiveMinusFixed),[],'omitnan');
    S.meanAdaptiveGridRatio = mean(T.AdaptiveGridRatio,'omitnan');
    S.medianAdaptiveGridRatio = median(T.AdaptiveGridRatio,'omitnan');
    S.meanFixedSeconds = mean(T.FixedSeconds,'omitnan');
    S.meanAdaptiveSeconds = mean(T.AdaptiveSeconds,'omitnan');
    S.medianRuntimeSpeedup = median(T.FixedSeconds./T.AdaptiveSeconds,'omitnan');
    S.maxValidatorSeconds = max(T.ValidatorSeconds,[],'omitnan');
    S.absErrorTolerance = tol;
    S.fixedAllPass = all(absF<=tol | isnan(absF));
    S.adaptiveAllPass = all(absA<=tol | isnan(absA));
    S.nRows = height(T);
    S.nAdaptiveFailures = sum(absA>tol);
    S.worstFixedRow = idxF;
    S.worstAdaptiveRow = idxA;
end


function print_summary(S,T)
    fprintf('\nPRODUCTION Y-GRID ACCURACY SUMMARY\n');
    fprintf('Rows evaluated                         : %d\n',S.nRows);
    fprintf('Tolerance                              : %.3e bit/symbol\n',S.absErrorTolerance);
    fprintf('Fixed mean/max |error|                 : %.3e / %.3e\n', ...
        S.meanAbsFixedError,S.maxAbsFixedError);
    fprintf('Adaptive mean/max |error|              : %.3e / %.3e\n', ...
        S.meanAbsAdaptiveError,S.maxAbsAdaptiveError);
    fprintf('Max |adaptive-fixed AMI|               : %.3e\n',S.maxAbsAdaptiveMinusFixed);
    fprintf('Adaptive all rows pass                 : %d\n',S.adaptiveAllPass);
    fprintf('Adaptive failures                      : %d\n',S.nAdaptiveFailures);
    fprintf('Median adaptive/fixed N_y ratio        : %.3f\n',S.medianAdaptiveGridRatio);
    fprintf('Median fixed/adaptive eval speedup     : %.3fx\n',S.medianRuntimeSpeedup);

    if ~isempty(S.worstAdaptiveRow) && isfinite(S.worstAdaptiveRow)
        r = T(S.worstAdaptiveRow,:);
        fprintf(['Worst adaptive case: M=%d SNR=%.1f sigma_X^2=%.3f geometry=%s ' ...
            '| xMax=%.4g | error=%+.3e | N_y=%d/%d\n'], ...
            r.M,r.SNRdB,r.SigmaX2,char(r.Geometry),r.XMax,r.AdaptiveError, ...
            r.AdaptiveYGridN,r.FixedYGridN);
    end
end


function validate_spacing(MVec,Pavg,dmin)
    req = (MVec-1)*dmin/2;
    bad = find(req>Pavg+100*eps(max(1,Pavg)),1,'first');
    if ~isempty(bad)
        error('diagnose_production_y_grid_accuracy:InfeasibleGap', ...
            'M=%d and d_min=%.6g require mean power %.6g > P_avg=%.6g.', ...
            MVec(bad),dmin,req(bad),Pavg);
    end
end


function tf = geometry_modes_valid(v)
    if isstring(v), v=cellstr(v); end
    tf = iscell(v) && ~isempty(v);
    if ~tf, return; end
    allowed = ["uniform","upper-quarter","upper-eighth","upper-one"];
    try
        q = lower(string(v));
        tf = all(ismember(q,allowed));
    catch
        tf = false;
    end
end

function t = sanitize_label(v)
    t = strtrim(char(string(v)));
    if isempty(t), t='run'; return; end
    t = regexprep(t,'[^A-Za-z0-9._-]+','_');
    t = regexprep(t,'_+','_');
    t = regexprep(t,'^[._-]+|[._-]+$','');
    if isempty(t), t='run'; end
end

function [d,id] = make_unique_dir(parent,id)
    if ~exist(parent,'dir'), mkdir(parent); end
    d = fullfile(parent,id); base=id; n=1;
    while exist(d,'dir')
        n=n+1; id=sprintf('%s_%02d',base,n); d=fullfile(parent,id);
    end
    [ok,msg]=mkdir(d);
    if ~ok, error('diagnose_production_y_grid_accuracy:mkdir','%s',msg); end
end

function token = num_token(v,nDec)
    token = sprintf(['%0.' num2str(nDec) 'f'],double(v));
    token = strrep(token,'-','m'); token = strrep(token,'+','p'); token = strrep(token,'.','p');
end

function tf=positive_scalar(v), tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&isreal(v)&&v>0; end
function tf=nonnegative_scalar(v), tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&isreal(v)&&v>=0; end
function tf=positive_integer(v), tf=positive_scalar(v)&&mod(v,1)==0; end
function tf=positive_integer_vector(v), tf=isnumeric(v)&&isvector(v)&&~isempty(v)&&all(isfinite(v))&&all(v>=1)&&all(mod(v,1)==0); end
function tf=finite_vector(v), tf=isnumeric(v)&&isvector(v)&&~isempty(v)&&all(isfinite(v)); end
function tf=nonnegative_vector(v), tf=finite_vector(v)&&all(v>=0); end
