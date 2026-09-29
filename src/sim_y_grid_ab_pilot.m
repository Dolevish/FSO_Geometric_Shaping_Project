function bundle = sim_y_grid_ab_pilot(varargin)
%SIM_Y_GRID_AB_PILOT  Commit-I.6 optimizer-level fixed/adaptive y-grid A/B.
%
% Controlled comparison of the production log-h fast AMI evaluator using:
%   A) historical fixed-bound y-grid;
%   B) candidate-adaptive y-grid with support scale 1.1 by default.
%
% The two methods use identical physical cases, SA hyperparameters, uniform
% PAM initialization, and replicate master seeds.  The adaptive grid changes
% numerical integration support only; it does not clip/project candidates or
% introduce a peak-power constraint.
%
% Scientific defaults:
%   M=[8 32], SNR=[20 30] dB, sigma_X^2=0.1, Pavg=1, d_min=0.01
%   log-h dt=0.01, span=+/-8 sigma_t, yBlock=512
%   adaptive scale=1.1
%   SA=4 starts x 3000 iterations, 2 independent replicates
%   T0=.4, Tf=1e-3, baseStd0=.15, itersPerTemp=50
%
% Suggested smoke test before the full pilot:
%   b = sim_y_grid_ab_pilot('MVec',4,'SNRVec',30, ...
%       'saMaxIter',100,'saNStarts',2,'nReplicates',1, ...
%       'UseParallelCases',false,'RunLabel','I6_smoke');

    p=inputParser;
    p.FunctionName=mfilename;
    addParameter(p,'MVec',[8 32],@positive_integer_vector);
    addParameter(p,'P_avg',1,@positive_scalar);
    addParameter(p,'SNRVec',[20 30],@(v) isnumeric(v)&&isvector(v)&&~isempty(v)&&all(isfinite(v)));
    addParameter(p,'TurbulenceVec',0.1,@nonnegative_vector);
    addParameter(p,'minGap',0.01,@nonnegative_scalar);
    addParameter(p,'pinZero',true,@(v) islogical(v)&&isscalar(v));

    addParameter(p,'loghDt',0.01,@positive_scalar);
    addParameter(p,'loghSpanSigma',8,@positive_scalar);
    addParameter(p,'loghYBlockSize',512,@positive_integer);
    addParameter(p,'AdaptiveScale',1.1,@(v) positive_scalar(v)&&v>=1);

    addParameter(p,'saMaxIter',3000,@positive_integer);
    addParameter(p,'saNStarts',4,@(v) positive_integer(v)&&v>=2);
    addParameter(p,'nReplicates',2,@positive_integer);
    addParameter(p,'T0',0.4,@positive_scalar);
    addParameter(p,'Tf',1e-3,@positive_scalar);
    addParameter(p,'BaseStd0',0.15,@positive_scalar);
    addParameter(p,'ItersPerTemp',50,@positive_integer);
    addParameter(p,'baseSeed',20260928,@nonnegative_scalar);
    addParameter(p,'historyEvery',50,@positive_integer);

    addParameter(p,'UseParallelCases',true,@(v) islogical(v)&&isscalar(v));
    addParameter(p,'RunLabel','',@(v) ischar(v)||(isstring(v)&&isscalar(v)));
    addParameter(p,'ResultsRoot',fso_result_utils.default_results_root(),@(v) ischar(v)||isstring(v));
    parse(p,varargin{:});
    o=p.Results;

    if ~(o.T0>o.Tf)
        error('sim_y_grid_ab_pilot:TemperatureOrder','T0 must exceed Tf.');
    end

    MVec=unique(double(o.MVec(:).'),'stable');
    snrVec=unique(double(o.SNRVec(:).'),'stable');
    sigVec=unique(double(o.TurbulenceVec(:).'),'stable');
    nRep=double(o.nReplicates);
    methods={'fixed-bound','candidate-adaptive'};
    nMethods=2;

    validate_spacing(MVec,double(o.P_avg),double(o.minGap));

    nM=numel(MVec); nSNR=numel(snrVec); nSig=numel(sigVec);
    [iMethod,iM,iSig,iSNR,iRep]=ndgrid(1:nMethods,1:nM,1:nSig,1:nSNR,1:nRep);
    tm=iMethod(:); tM=iM(:); ts=iSig(:); tn=iSNR(:); tr=iRep(:);
    nTasks=numel(tm);

    resultsRoot=char(o.ResultsRoot);
    parentDir=fullfile(resultsRoot,'y_grid_ab_pilot');
    label=sanitize_label(o.RunLabel);
    configToken=sprintf('LOGHdt%s_T%s_scale%s_I%d_S%d_R%d', ...
        num_token(o.loghDt,4),num_token(o.loghSpanSigma,2), ...
        num_token(o.AdaptiveScale,2),o.saMaxIter,o.saNStarts,nRep);
    utcToken=char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmssSSS''Z'''));
    requested=sprintf('%s_%s_%s',label,configToken,utcToken);
    [outDir,runId]=make_unique_dir(parentDir,requested);
    caseDir=fullfile(outDir,'cases'); mkdir(caseDir);

    runMeta=fso_result_utils.run_metadata(mfilename);
    runMeta.runId=runId;
    runMeta.runLabel=char(string(o.RunLabel));
    runMeta.outputDirectory=outDir;

    fprintf('\n============================================================\n');
    fprintf('Commit I.6 - fixed-bound vs candidate-adaptive y-grid SA A/B\n');
    fprintf('Run ID: %s\n',runId);
    fprintf('M=%s | SNR=%s dB | sigma_X^2=%s | reps=%d\n', ...
        mat2str(MVec),mat2str(snrVec),mat2str(sigVec),nRep);
    fprintf('log-h dt=%.5g | span=+/-%.3g sigma_t | yBlock=%d\n', ...
        o.loghDt,o.loghSpanSigma,o.loghYBlockSize);
    fprintf('Adaptive scale=%.4g | fixed mode retained as control\n',o.AdaptiveScale);
    fprintf('SA=%d starts x %d iterations | T0=%.4g Tf=%.4g step0=%.4g block=%d\n', ...
        o.saNStarts,o.saMaxIter,o.T0,o.Tf,o.BaseStd0,o.ItersPerTemp);
    fprintf('Tasks=%d (2 methods x physical cases x replicates)\n',nTasks);
    fprintf('Results: %s\n',outDir);
    fprintf('============================================================\n\n');

    caseResults=cell(nTasks,1); caseFiles=cell(nTasks,1);
    doPar=logical(o.UseParallelCases);
    if doPar
        try, doPar=license('test','Distrib_Computing_Toolbox')~=0; catch, doPar=false; end
    end

    if doPar
        pool=gcp('nocreate'); if isempty(pool), pool=parpool('Processes'); end
        fprintf('Outer case parallelism: %d process workers\n',pool.NumWorkers);
        Q=parallel.pool.DataQueue; doneCount=0; afterEach(Q,@(~)tick());
        parfor k=1:nTasks
            [caseResults{k},caseFiles{k}]=run_one(methods{tm(k)}, ...
                MVec(tM(k)),snrVec(tn(k)),sigVec(ts(k)),tr(k),o,runMeta,caseDir);
            send(Q,1);
        end
    else
        fprintf('Outer case parallelism: disabled\n');
        for k=1:nTasks
            [caseResults{k},caseFiles{k}]=run_one(methods{tm(k)}, ...
                MVec(tM(k)),snrVec(tn(k)),sigVec(ts(k)),tr(k),o,runMeta,caseDir);
            fprintf('  completed %d/%d\n',k,nTasks);
        end
    end

    summary=build_summary(caseResults,methods,MVec,sigVec,snrVec,nRep);

    bundle=struct();
    bundle.meta=runMeta;
    bundle.experiment='y_grid_ab_pilot';
    bundle.run=struct('id',runId,'label',char(string(o.RunLabel)), ...
        'outputDirectory',outDir);
    bundle.axes=struct('method',{methods},'M',MVec,'P_avg',double(o.P_avg), ...
        'SNR_dB',snrVec,'sigma_X_sq',sigVec,'replicate',1:nRep);
    bundle.method=struct('fastFadingMethod','logh','loghDt',double(o.loghDt), ...
        'loghSpanSigma',double(o.loghSpanSigma),'loghYBlockSize',double(o.loghYBlockSize), ...
        'fixedYGridMode','fixed-bound','adaptiveYGridMode','candidate-adaptive', ...
        'adaptiveScale',double(o.AdaptiveScale));
    bundle.sa=struct('maxIter',double(o.saMaxIter),'nStarts',double(o.saNStarts), ...
        'nReplicates',nRep,'T0',double(o.T0),'Tf',double(o.Tf), ...
        'baseStd0',double(o.BaseStd0),'itersPerTemp',double(o.ItersPerTemp), ...
        'baseSeed',double(o.baseSeed),'historyEvery',double(o.historyEvery));
    bundle.caseFiles=caseFiles;
    bundle.summary=summary;

    bundlePath=fullfile(outDir,'yGridABPilot.mat');
    save(bundlePath,'bundle','-v7.3');
    fprintf('\nSaved I.6 A/B bundle: %s\n',bundlePath);
    print_summary(summary);

    function tick()
        doneCount=doneCount+1;
        fprintf('  completed %d/%d\n',doneCount,nTasks);
    end
end

function [c,path]=run_one(mode,M,SNR_dB,sigma_X_sq,rep,o,runMeta,caseDir)
    seedBase=double(o.baseSeed)+(double(rep)-1)*10000019;
    seed=fso_result_utils.case_seed(seedBase,M,SNR_dB,sigma_X_sq);
    xMaxBound=fso_result_utils.feasible_xmax_bound(M,o.P_avg,o.minGap,o.pinZero);

    cfg=build_fso_config(M,o.P_avg,SNR_dB,sigma_X_sq, ...
        'FastFadingMethod','logh','loghDt',o.loghDt, ...
        'loghSpanSigma',o.loghSpanSigma,'loghYBlockSize',o.loghYBlockSize, ...
        'YGridMode',mode,'YGridScale',o.AdaptiveScale,'xMaxBound',xMaxBound, ...
        'saMaxIter',o.saMaxIter,'saNStarts',o.saNStarts,'saUseParallel',false, ...
        'minGap',o.minGap,'pinZero',o.pinZero,'seedInit',seed, ...
        'logEvery',0,'historyEvery',o.historyEvery);

    cfg.SA.T0=double(o.T0); cfg.SA.Tf=double(o.Tf);
    cfg.SA.baseStd0=double(o.BaseStd0);
    cfg.SA.itersPerTemp=double(o.ItersPerTemp);
    cfg.SA.nBlocks=ceil(cfg.SA.maxIter/cfg.SA.itersPerTemp);
    cfg.SA.coolingRate=exp(log(cfg.SA.Tf/cfg.SA.T0)/cfg.SA.nBlocks);

    xPAM=define_constellation(M,o.P_avg,true,"mean");
    AMI_functions.assert_constellation_feasible(xPAM,cfg,1e-10);
    baselineVal=cfg.AMI_Validator(xPAM(:));
    baselineFast=cfg.AMI_Evaluator(xPAM(:));

    label=sprintf('[I6|%s|M%d|SNR%.1f|sig%.3f|rep%d]', ...
        short_mode(mode),M,SNR_dB,sigma_X_sq,rep);
    t=tic;
    [saOut,starts]=sa_multistart(cfg,xPAM,label);
    wall=toc(t);

    win=starts(saOut.bestStartValidated);
    restartVals=[starts.bestMIValidated];
    k99=convergence_iteration(win.history,win.initialMI,win.bestMIFast,0.99);
    nFast=double(o.saNStarts)*double(o.saMaxIter);

    c=struct();
    c.meta=runMeta;
    c.case=struct('method',mode,'M',M,'P_avg',double(o.P_avg),'SNR_dB',SNR_dB, ...
        'sigma_X_sq',sigma_X_sq,'replicate',rep,'seed',seed);
    c.config=fso_result_utils.config_snapshot(cfg);
    c.baseline=struct('x',xPAM(:),'amiValidated',baselineVal,'amiFast',baselineFast, ...
        'fastValidationGap',baselineFast-baselineVal);
    c.winner=struct('x',saOut.bestXValidated(:),'amiValidated',saOut.bestMIValidated, ...
        'amiFast',saOut.validatedWinnerFastMI,'gainBits',saOut.bestMIValidated-baselineVal, ...
        'k99',k99,'acceptanceRate',win.acceptanceRate);
    c.validation=struct('selectionChanged',logical(saOut.selectionChanged), ...
        'meanAbsFastValidationGap',saOut.meanAbsFastValidationGap, ...
        'maxAbsFastValidationGap',saOut.maxAbsFastValidationGap);
    c.repeatability=struct('restartValidatedMean',mean(restartVals,'omitnan'), ...
        'restartValidatedStd',std(restartVals,0,'omitnan'), ...
        'restartValidatedRange',max(restartVals)-min(restartVals));
    c.runtime=struct('wallSeconds',wall,'saSecondsSummed',saOut.totalSARuntime, ...
        'validationSeconds',saOut.validationTotalRuntime,'fastEvaluations',nFast, ...
        'fastEvaluationsPerSecond',nFast/max(saOut.totalSARuntime,eps));
    c.geometry=struct('xMax',max(saOut.bestXValidated), ...
        'xMaxBound',xMaxBound,'xMaxToBoundRatio',max(saOut.bestXValidated)/xMaxBound);
    c.starts=starts;

    fname=sprintf('case_%s_M%d_SNR%s_sig%s_rep%02d.mat', ...
        short_mode(mode),M,num_token(SNR_dB,2),num_token(sigma_X_sq,4),rep);
    path=fullfile(caseDir,fname);
    save(path,'c','-v7.3');
end

function S=build_summary(cases,methods,MVec,sigVec,snrVec,nRep)
    rows=cell(numel(cases),16);
    for k=1:numel(cases)
        c=cases{k};
        rows(k,:)={c.case.method,c.case.M,c.case.sigma_X_sq,c.case.SNR_dB,c.case.replicate, ...
            c.winner.amiValidated,c.winner.amiFast,c.winner.gainBits,c.winner.k99, ...
            c.validation.maxAbsFastValidationGap,double(c.validation.selectionChanged), ...
            c.repeatability.restartValidatedStd,c.runtime.wallSeconds, ...
            c.runtime.fastEvaluationsPerSecond,c.geometry.xMax,c.geometry.xMaxToBoundRatio};
    end
    T=cell2table(rows,'VariableNames',{'Method','M','SigmaX2','SNRdB','Replicate', ...
        'AMIValidated','AMIFast','GainBits','K99','MaxAbsFastValidationGap', ...
        'SelectionChanged','RestartValidatedStd','WallSeconds','FastEvalPerSecond', ...
        'XMax','XMaxToBoundRatio'});
    T.Method=string(T.Method);
    numericNames=setdiff(T.Properties.VariableNames,{'Method'});
    for j=1:numel(numericNames)
        name=numericNames{j};
        if iscell(T.(name))
            T.(name)=cell2mat(T.(name));
        end
    end

    pairRows={}; q=0;
    for iM=1:numel(MVec)
        for iSig=1:numel(sigVec)
            for iS=1:numel(snrVec)
                for r=1:nRep
                    key=T.M==MVec(iM)&T.SigmaX2==sigVec(iSig)&T.SNRdB==snrVec(iS)&T.Replicate==r;
                    f=T(key&T.Method=="fixed-bound",:);
                    a=T(key&T.Method=="candidate-adaptive",:);
                    if height(f)~=1||height(a)~=1, error('sim_y_grid_ab_pilot:Pairing','Missing A/B pair.'); end
                    q=q+1;
                    pairRows(q,:)={MVec(iM),sigVec(iSig),snrVec(iS),r, ...
                        f.AMIValidated,a.AMIValidated,a.AMIValidated-f.AMIValidated, ...
                        f.WallSeconds,a.WallSeconds,f.WallSeconds/a.WallSeconds, ...
                        f.K99,a.K99,f.RestartValidatedStd,a.RestartValidatedStd, ...
                        f.MaxAbsFastValidationGap,a.MaxAbsFastValidationGap}; %#ok<AGROW>
                end
            end
        end
    end
    P=cell2table(pairRows,'VariableNames',{'M','SigmaX2','SNRdB','Replicate', ...
        'FixedAMI','AdaptiveAMI','AdaptiveMinusFixedAMI','FixedWallSeconds', ...
        'AdaptiveWallSeconds','WallSpeedupFixedOverAdaptive','FixedK99','AdaptiveK99', ...
        'FixedRestartStd','AdaptiveRestartStd','FixedMaxFVGap','AdaptiveMaxFVGap'});
    for j=1:numel(P.Properties.VariableNames)
        name=P.Properties.VariableNames{j};
        if iscell(P.(name))
            P.(name)=cell2mat(P.(name));
        end
    end

    S=struct(); S.raw=T; S.paired=P;
    S.overall=struct();
    S.overall.medianWallSpeedupFixedOverAdaptive=median(P.WallSpeedupFixedOverAdaptive,'omitnan');
    S.overall.meanWallSpeedupFixedOverAdaptive=mean(P.WallSpeedupFixedOverAdaptive,'omitnan');
    S.overall.maxAbsPairedValidatedAMIDifference=max(abs(P.AdaptiveMinusFixedAMI),[],'omitnan');
    S.overall.meanAdaptiveMinusFixedAMI=mean(P.AdaptiveMinusFixedAMI,'omitnan');
    S.overall.maxAdaptiveFastValidationGap=max(P.AdaptiveMaxFVGap,[],'omitnan');
    S.overall.adaptiveSelectionChangedFraction=mean(T.SelectionChanged(T.Method=="candidate-adaptive"),'omitnan');
    S.overall.fixedSelectionChangedFraction=mean(T.SelectionChanged(T.Method=="fixed-bound"),'omitnan');

    aggRows={}; q=0;
    for m=1:numel(methods)
        for iM=1:numel(MVec)
            for iSig=1:numel(sigVec)
                for iS=1:numel(snrVec)
                    idx=T.Method==string(methods{m})&T.M==MVec(iM)&T.SigmaX2==sigVec(iSig)&T.SNRdB==snrVec(iS);
                    q=q+1;
                    aggRows(q,:)={methods{m},MVec(iM),sigVec(iSig),snrVec(iS), ...
                        mean(T.AMIValidated(idx),'omitnan'),std(T.AMIValidated(idx),0,'omitnan'), ...
                        mean(T.GainBits(idx),'omitnan'),mean(T.K99(idx),'omitnan'), ...
                        mean(T.RestartValidatedStd(idx),'omitnan'), ...
                        max(T.MaxAbsFastValidationGap(idx),[],'omitnan'), ...
                        mean(T.SelectionChanged(idx),'omitnan'),mean(T.WallSeconds(idx),'omitnan'), ...
                        mean(T.FastEvalPerSecond(idx),'omitnan')}; %#ok<AGROW>
                end
            end
        end
    end
    A=cell2table(aggRows,'VariableNames',{'Method','M','SigmaX2','SNRdB', ...
        'MeanAMIValidated','StdAMIValidated','MeanGainBits','MeanK99', ...
        'MeanRestartValidatedStd','MaxFastValidationGap','SelectionChangedFraction', ...
        'MeanWallSeconds','MeanFastEvalPerSecond'});
    A.Method=string(A.Method);
    for j=2:numel(A.Properties.VariableNames)
        name=A.Properties.VariableNames{j};
        if iscell(A.(name))
            A.(name)=cell2mat(A.(name));
        end
    end
    S.aggregate=A;
end

function print_summary(S)
    fprintf('\nI.6 A/B SUMMARY\n');
    fprintf('median fixed/adaptive wall-time speedup: %.3fx\n',S.overall.medianWallSpeedupFixedOverAdaptive);
    fprintf('mean fixed/adaptive wall-time speedup  : %.3fx\n',S.overall.meanWallSpeedupFixedOverAdaptive);
    fprintf('max |paired validated AMI difference|  : %.3e bit/symbol\n',S.overall.maxAbsPairedValidatedAMIDifference);
    fprintf('mean adaptive-fixed validated AMI      : %+.3e bit/symbol\n',S.overall.meanAdaptiveMinusFixedAMI);
    fprintf('max adaptive fast-validation gap       : %.3e bit/symbol\n',S.overall.maxAdaptiveFastValidationGap);
    fprintf('selection changed fraction fixed/adapt : %.3f / %.3f\n', ...
        S.overall.fixedSelectionChangedFraction,S.overall.adaptiveSelectionChangedFraction);

    P=S.paired;
    fprintf('\nPer paired case:\n');
    for k=1:height(P)
        fprintf(['  M=%d SNR=%.1f sig=%.3f rep=%d | AMI fixed/adapt %.6f/%.6f ' ...
            '| dAMI=%+.2e | speedup %.2fx | k99 %.0f/%.0f\n'], ...
            P.M(k),P.SNRdB(k),P.SigmaX2(k),P.Replicate(k),P.FixedAMI(k),P.AdaptiveAMI(k), ...
            P.AdaptiveMinusFixedAMI(k),P.WallSpeedupFixedOverAdaptive(k),P.FixedK99(k),P.AdaptiveK99(k));
    end
end

function k=convergence_iteration(history,initialMI,bestMI,fraction)
    improvement=bestMI-initialMI;
    if ~(isfinite(improvement)&&improvement>0), k=0; return; end
    target=initialMI+fraction*improvement;
    idx=find(history.bestMI>=target,1,'first');
    if isempty(idx), k=history.iter(end); else, k=history.iter(idx); end
end

function validate_spacing(MVec,Pavg,d)
    req=(MVec-1)*d/2;
    bad=find(req>Pavg+100*eps(max(1,Pavg)),1);
    if ~isempty(bad)
        error('sim_y_grid_ab_pilot:InfeasibleGap', ...
            'M=%d with d_min=%.6g requires mean power %.6g > %.6g.',MVec(bad),d,req(bad),Pavg);
    end
end

function s=short_mode(mode)
    if strcmp(mode,'fixed-bound'), s='fixed'; else, s='adaptive'; end
end

function t=sanitize_label(v)
    t=strtrim(char(string(v))); if isempty(t),t='run';return;end
    t=regexprep(t,'[^A-Za-z0-9._-]+','_'); t=regexprep(t,'_+','_');
    t=regexprep(t,'^[._-]+|[._-]+$',''); if isempty(t),t='run';end
end

function [d,id]=make_unique_dir(parent,id)
    if ~exist(parent,'dir'),mkdir(parent);end
    d=fullfile(parent,id); n=1; base=id;
    while exist(d,'dir'), n=n+1; id=sprintf('%s_%02d',base,n); d=fullfile(parent,id); end
    [ok,msg]=mkdir(d); if ~ok,error('sim_y_grid_ab_pilot:mkdir','%s',msg);end
end

function token=num_token(v,nDec)
    token=sprintf(['%0.' num2str(nDec) 'f'],double(v));
    token=strrep(token,'-','m'); token=strrep(token,'+','p'); token=strrep(token,'.','p');
end

function tf=positive_scalar(v), tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&isreal(v)&&v>0; end
function tf=nonnegative_scalar(v), tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&isreal(v)&&v>=0; end
function tf=positive_integer(v), tf=positive_scalar(v)&&mod(v,1)==0; end
function tf=positive_integer_vector(v), tf=isnumeric(v)&&isvector(v)&&~isempty(v)&&all(isfinite(v))&&all(v>=1)&&all(mod(v,1)==0); end
function tf=nonnegative_vector(v), tf=isnumeric(v)&&isvector(v)&&~isempty(v)&&all(isfinite(v))&&all(v>=0); end