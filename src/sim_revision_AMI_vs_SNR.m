function bundle = sim_revision_AMI_vs_SNR(varargin)
%SIM_REVISION_AMI_VS_SNR  Production AMI-vs-SNR pipeline for IEEE revision.
%
% Commit I.9 production integration:
%   - uniform quadrature in t=ln(h)
%   - candidate-adaptive outer y-grid, scale=1.1
%   - deterministic per-channel-case log-h spacing policy
%   - independent validator and validated restart selection unchanged
%
% I.8 showed that dt=0.01 is sufficient over most of the production grid,
% but high-SNR/high-order stress cases require dt=0.005. I.9 therefore
% selects dt ONCE from (M,SNR,sigma_X^2) before a case starts and keeps it
% fixed throughout all SA iterations, restarts, and replicates of that
% physical channel point. The selection never depends on the current
% candidate constellation.
%
% Default I.9 dt policy:
%   base dt    = 0.01
%   refined dt = 0.005
%   refine if SNR>=30 dB and
%      ((M>=32 && sigma_X^2>=0.1) || (M>=16 && sigma_X^2>=0.3)).
% On the current production grid this selects exactly the four I.8 stress
% cases: (32,30,.1), (32,30,.2), (32,30,.3), (16,30,.3).
%
% Backward compatibility: if 'loghDt' is supplied explicitly while
% 'LoghDtPolicy' is omitted, the run is interpreted as a fixed-dt historical
% run. To force a fixed value explicitly use LoghDtPolicy='fixed'.
%
% Scientific defaults:
%   M                  = [4 8 16 32]
%   P_avg              = 1
%   SNR                 = 5:5:30 dB
%   sigma_X^2           = [0 0.1 0.2 0.3]
%   d_min               = 0.01
%   fast fading method  = uniform log-h quadrature
%   log-h dt policy     = case-adaptive (0.01 / 0.005)
%   log-h span          = mu_t +/- 8 sigma_t
%   y-grid mode         = candidate-adaptive
%   y-grid scale        = 1.1
%
% Name-value options include:
%   'LoghDtPolicy'      'case-adaptive' (default) or 'fixed'
%   'loghDt'            base dt, or fixed dt in fixed mode (default 0.01)
%   'loghDtRefined'     refined dt in case-adaptive mode (default 0.005)
% plus the historical options from Commit I.7.

    p = inputParser;
    p.FunctionName = mfilename;

    addParameter(p,'MVec',[4 8 16 32],@(v) isnumeric(v) && isvector(v) && ~isempty(v) && all(isfinite(v)) && all(v>=2) && all(mod(v,1)==0));
    addParameter(p,'P_avg',1,@positive_scalar);
    addParameter(p,'SNRVec',5:5:30,@(v) isnumeric(v) && isvector(v) && ~isempty(v) && all(isfinite(v)));
    addParameter(p,'TurbulenceVec',[0 0.1 0.2 0.3],@(v) isnumeric(v) && isvector(v) && ~isempty(v) && all(isfinite(v)) && all(v>=0));
    addParameter(p,'minGap',0.01,@(v) isnumeric(v) && isscalar(v) && isfinite(v) && isreal(v) && v>=0);
    addParameter(p,'pinZero',true,@(v) islogical(v) && isscalar(v));
    addParameter(p,'FastFadingMethod','logh',@(v) ischar(v) || (isstring(v) && isscalar(v)));
    addParameter(p,'ghN_h',400,@positive_integer);
    addParameter(p,'LoghDtPolicy','case-adaptive',@(v) ischar(v) || (isstring(v) && isscalar(v)));
    addParameter(p,'loghDt',0.01,@positive_scalar);
    addParameter(p,'loghDtRefined',0.005,@positive_scalar);
    addParameter(p,'loghSpanSigma',8,@positive_scalar);
    addParameter(p,'loghYBlockSize',512,@positive_integer);
    addParameter(p,'YGridMode','candidate-adaptive',@(v) ischar(v) || (isstring(v) && isscalar(v)));
    addParameter(p,'YGridScale',1.1,@(v) isnumeric(v) && isscalar(v) && isfinite(v) && isreal(v) && v>=1);

    addParameter(p,'saMaxIter',8000,@positive_integer);
    addParameter(p,'saNStarts',6,@positive_integer);
    addParameter(p,'nReplicates',3,@positive_integer);
    addParameter(p,'T0',0.4,@positive_scalar);
    addParameter(p,'Tf',1e-3,@positive_scalar);
    addParameter(p,'BaseStd0',0.15,@positive_scalar);
    addParameter(p,'ItersPerTemp',50,@positive_integer);
    addParameter(p,'baseSeed',20260928,@(v) isnumeric(v) && isscalar(v) && isfinite(v) && v>=0);
    addParameter(p,'historyEvery',50,@positive_integer);

    addParameter(p,'UseParallelCases',true,@(v) islogical(v) && isscalar(v));
    addParameter(p,'MakePlots',true,@(v) islogical(v) && isscalar(v));
    addParameter(p,'SaveFigures',true,@(v) islogical(v) && isscalar(v));
    addParameter(p,'RunLabel','',@(v) ischar(v) || (isstring(v) && isscalar(v)));
    addParameter(p,'ResultsRoot',fso_result_utils.default_results_root(),@(v) ischar(v) || isstring(v));
    parse(p,varargin{:});
    o = p.Results;

    policyExplicit = ~ismember('LoghDtPolicy',p.UsingDefaults);
    dtExplicit = ~ismember('loghDt',p.UsingDefaults);
    if dtExplicit && ~policyExplicit
        % Preserve historical scripts that explicitly supplied loghDt.
        o.LoghDtPolicy = 'fixed';
    end

    fastMethod = lower(string(o.FastFadingMethod));
    if ~ismember(fastMethod,["logh","gh"])
        error('sim_revision_AMI_vs_SNR:BadFastFadingMethod', ...
            'FastFadingMethod must be ''logh'' or ''gh'', got ''%s''.',fastMethod);
    end
    yGridMode = lower(string(o.YGridMode));
    if ~ismember(yGridMode,["fixed-bound","candidate-adaptive"])
        error('sim_revision_AMI_vs_SNR:BadYGridMode', ...
            'YGridMode must be ''fixed-bound'' or ''candidate-adaptive'', got ''%s''.',yGridMode);
    end
    dtPolicy = lower(string(o.LoghDtPolicy));
    if ~ismember(dtPolicy,["case-adaptive","fixed"])
        error('sim_revision_AMI_vs_SNR:BadLoghDtPolicy', ...
            'LoghDtPolicy must be ''case-adaptive'' or ''fixed'', got ''%s''.',dtPolicy);
    end
    if o.loghDtRefined > o.loghDt
        error('sim_revision_AMI_vs_SNR:BadLoghDtOrder', ...
            'loghDtRefined must not exceed loghDt.');
    end
    if fastMethod=="gh" && yGridMode=="candidate-adaptive"
        error('sim_revision_AMI_vs_SNR:AdaptiveYGridRequiresLogh', ...
            ['Legacy GH production diagnostics require YGridMode=''fixed-bound''. ' ...
             'candidate-adaptive is validated only for the log-h evaluator.']);
    end
    o.FastFadingMethod = char(fastMethod);
    o.YGridMode = char(yGridMode);
    o.LoghDtPolicy = char(dtPolicy);

    if ~(double(o.T0) > double(o.Tf))
        error('sim_revision_AMI_vs_SNR:TemperatureOrder','T0 must be greater than Tf.');
    end

    MVec = unique(double(o.MVec(:).'),'stable');
    snrVec = unique(double(o.SNRVec(:).'),'stable');
    sigVec = unique(double(o.TurbulenceVec(:).'),'stable');
    P_avg = double(o.P_avg);
    nRep = double(o.nReplicates);
    validate_spacing_feasibility(MVec,P_avg,double(o.minGap));

    nM = numel(MVec); nSig = numel(sigVec); nSNR = numel(snrVec);
    nPhysical = nM*nSig*nSNR; nCases = nPhysical*nRep;

    xMaxBounds = nan(1,nM);
    for iM = 1:nM
        xMaxBounds(iM) = fso_result_utils.feasible_xmax_bound( ...
            MVec(iM),P_avg,double(o.minGap),logical(o.pinZero));
    end

    % Commit I.9: resolve the integration spacing once per physical case.
    dtGrid = repmat(double(o.loghDt),nM,nSig,nSNR);
    dtMeta = cell(nM,nSig,nSNR);
    if fastMethod=="logh"
        for iM=1:nM
            for iSig=1:nSig
                for iSNR=1:nSNR
                    [dtGrid(iM,iSig,iSNR),dtMeta{iM,iSig,iSNR}] = select_logh_dt( ...
                        MVec(iM),snrVec(iSNR),sigVec(iSig), ...
                        'Policy',dtPolicy,'BaseDt',o.loghDt, ...
                        'RefinedDt',o.loghDtRefined,'FixedDt',o.loghDt);
                end
            end
        end
    else
        for k=1:numel(dtMeta)
            dtMeta{k}=struct('policy','not-applicable','ruleVersion','I9-v1', ...
                'selectedDt',NaN,'baseDt',NaN,'refinedDt',NaN,'fixedDt',NaN, ...
                'isRefined',false,'hardCaseRuleTriggered',false, ...
                'reason','Gauss-Hermite evaluator: log-h dt not applicable');
        end
    end
    refinedMask = cellfun(@(m) isfield(m,'isRefined') && m.isRefined,dtMeta);
    nRefinedPhysical = sum(refinedMask(:));

    resultsRoot = char(o.ResultsRoot);
    if fastMethod == "logh"
        if yGridMode=="candidate-adaptive", gridTag=sprintf('Yadapt%s',num_token(o.YGridScale,2));
        else, gridTag='Yfixed'; end
        if dtPolicy=="case-adaptive"
            dtTag=sprintf('LOGHdtAUTO%s-%s',num_token(o.loghDt,4),num_token(o.loghDtRefined,4));
        else
            dtTag=sprintf('LOGHdt%s',num_token(o.loghDt,4));
        end
        methodTag = sprintf('gap%s_%s_T%s_%s',num_token(o.minGap,4),dtTag, ...
            num_token(o.loghSpanSigma,2),gridTag);
    else
        methodTag = sprintf('gap%s_GH%d',num_token(o.minGap,4),o.ghN_h);
    end
    parentDir = fullfile(resultsRoot,'revised_ami_vs_snr',methodTag);

    runMeta = fso_result_utils.run_metadata(mfilename);
    labelToken = sanitize_run_label(o.RunLabel);
    configToken = sprintf('I%d_S%d_R%d',o.saMaxIter,o.saNStarts,nRep);
    utcToken = char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmssSSS''Z'''));
    requestedRunId = sprintf('%s_%s_%s',labelToken,configToken,utcToken);
    [outDir,runId] = make_unique_run_directory(parentDir,requestedRunId);
    caseDir = fullfile(outDir,'cases'); figDir = fullfile(outDir,'figures');
    mkdir(caseDir); if o.MakePlots && o.SaveFigures, mkdir(figDir); end

    runMeta.runId = runId; runMeta.runLabel = char(string(o.RunLabel));
    runMeta.outputDirectory = outDir; runMeta.methodTag = methodTag;

    saProfile = struct('T0',double(o.T0),'Tf',double(o.Tf), ...
        'baseStd0',double(o.BaseStd0),'itersPerTemp',double(o.ItersPerTemp));

    fprintf('\n============================================================\n');
    fprintf('IEEE revision production AMI-vs-SNR pipeline (Commit I.9)\n');
    fprintf('Run ID: %s\n',runId);
    if strlength(string(o.RunLabel))>0, fprintf('Run label: %s\n',char(string(o.RunLabel))); end
    fprintf('M=%s | SNR[dB]=%s | sigma_X^2=%s | Pavg=%.4g\n', ...
        mat2str(MVec),mat2str(snrVec),mat2str(sigVec),P_avg);
    fprintf('Revised geometry: d_min=%.6g | pinZero=%d\n',o.minGap,o.pinZero);
    if fastMethod=="logh"
        if dtPolicy=="case-adaptive"
            fprintf('Fast fading integration: uniform log-h | dt policy=case-adaptive | base=%.6g | refined=%.6g\n', ...
                o.loghDt,o.loghDtRefined);
            fprintf('I.9 refined physical cases: %d/%d\n',nRefinedPhysical,nPhysical);
        else
            fprintf('Fast fading integration: uniform log-h | fixed dt=%.6g\n',o.loghDt);
        end
        fprintf('log-h span=+/-%.3g sigma_t | yBlock=%d\n',o.loghSpanSigma,o.loghYBlockSize);
        fprintf('Outer y-grid: %s | adaptive scale=%.4g\n',char(yGridMode),o.YGridScale);
    else
        fprintf('Fast fading integration: Gauss-Hermite | order=%d | y-grid=fixed-bound\n',o.ghN_h);
    end
    fprintf('SA=%d restarts x %d iterations | independent replicates=%d\n',o.saNStarts,o.saMaxIter,nRep);
    fprintf('SA profile: T0=%.5g | Tf=%.5g | step0=%.5g | block=%d\n', ...
        saProfile.T0,saProfile.Tf,saProfile.baseStd0,saProfile.itersPerTemp);
    fprintf('Physical channel points=%d | optimization cases=%d\n',nPhysical,nCases);
    fprintf('Validated restart selection: enabled\n');
    fprintf('Results: %s\n',outDir);
    fprintf('============================================================\n\n');

    baseline = build_baseline_grid(MVec,sigVec,snrVec,P_avg,o,xMaxBounds,dtGrid,dtMeta);

    [iMGrid,iSigGrid,iSNRGrid,iRGrid] = ndgrid(1:nM,1:nSig,1:nSNR,1:nRep);
    taskM=iMGrid(:); taskSig=iSigGrid(:); taskSNR=iSNRGrid(:); taskR=iRGrid(:);
    caseResults=cell(nCases,1); caseFiles=cell(nCases,1);

    doPar=logical(o.UseParallelCases);
    if doPar
        try, doPar=license('test','Distrib_Computing_Toolbox')~=0; catch, doPar=false; end
    end
    if doPar
        pool=gcp('nocreate'); if isempty(pool), pool=parpool('Processes'); end
        fprintf('Outer case parallelism: %d process workers\n',pool.NumWorkers);
        Q=parallel.pool.DataQueue; doneCount=0; afterEach(Q,@(~) progress_tick());
        parfor k=1:nCases
            b=baseline{taskM(k),taskSig(k),taskSNR(k)};
            [caseResults{k},caseFiles{k}] = run_revision_case( ...
                taskM(k),taskSig(k),taskSNR(k),taskR(k),MVec,sigVec,snrVec, ...
                P_avg,o,xMaxBounds,saProfile,b,runMeta,caseDir,dtGrid,dtMeta);
            send(Q,1);
        end
    else
        fprintf('Outer case parallelism: disabled\n');
        for k=1:nCases
            b=baseline{taskM(k),taskSig(k),taskSNR(k)};
            [caseResults{k},caseFiles{k}] = run_revision_case( ...
                taskM(k),taskSig(k),taskSNR(k),taskR(k),MVec,sigVec,snrVec, ...
                P_avg,o,xMaxBounds,saProfile,b,runMeta,caseDir,dtGrid,dtMeta);
            fprintf('  completed %d/%d\n',k,nCases);
        end
    end

    caseGrid=reshape(caseResults,nM,nSig,nSNR,nRep);
    fileGrid=reshape(caseFiles,nM,nSig,nSNR,nRep);
    summary=build_summary(caseGrid,baseline,MVec,sigVec,snrVec,nRep);

    bundle=struct();
    bundle.meta=runMeta; bundle.experiment='revised_AMI_vs_SNR';
    bundle.run=struct('id',runId,'label',char(string(o.RunLabel)), ...
        'methodTag',methodTag,'outputDirectory',outDir);
    bundle.axes=struct('M',MVec,'P_avg',P_avg,'SNR_dB',snrVec, ...
        'sigma_X_sq',sigVec,'replicate',1:nRep);
    bundle.method=struct('minGap',double(o.minGap),'pinZero',logical(o.pinZero), ...
        'fastFadingMethod',char(fastMethod),'ghN_h',double(o.ghN_h), ...
        'loghDt',double(o.loghDt),'loghDtPolicy',char(dtPolicy), ...
        'loghDtBase',double(o.loghDt),'loghDtRefined',double(o.loghDtRefined), ...
        'loghDtRuleVersion','I9-v1','loghDtGrid',dtGrid, ...
        'nRefinedPhysicalCases',nRefinedPhysical, ...
        'loghSpanSigma',double(o.loghSpanSigma),'loghYBlockSize',double(o.loghYBlockSize), ...
        'yGridMode',char(yGridMode),'yGridScale',double(o.YGridScale), ...
        'validatedRestartSelection',true,'independentReplicates',true,'outerParallelism',doPar);
    bundle.sa=struct('maxIter',double(o.saMaxIter),'nStarts',double(o.saNStarts), ...
        'nReplicates',nRep,'T0',saProfile.T0,'Tf',saProfile.Tf, ...
        'baseStd0',saProfile.baseStd0,'itersPerTemp',saProfile.itersPerTemp, ...
        'baseSeed',double(o.baseSeed),'historyEvery',double(o.historyEvery));
    bundle.xMaxBounds=xMaxBounds; bundle.baseline=baseline; bundle.caseFiles=fileGrid;
    bundle.summary=summary; bundle.figureFiles={};

    if o.MakePlots, bundle.figureFiles=make_summary_plots(bundle,figDir,logical(o.SaveFigures)); end
    bundlePath=fullfile(outDir,'revisedAMIvsSNR.mat'); save(bundlePath,'bundle','-v7.3');
    fprintf('\nSaved revised AMI-vs-SNR bundle: %s\n',bundlePath); print_summary(bundle);

    function progress_tick()
        doneCount=doneCount+1; fprintf('  completed %d/%d\n',doneCount,nCases);
    end
end


function baseline=build_baseline_grid(MVec,sigVec,snrVec,P_avg,o,xMaxBounds,dtGrid,dtMeta)
    nM=numel(MVec); nSig=numel(sigVec); nSNR=numel(snrVec);
    baseline=cell(nM,nSig,nSNR);
    fprintf('Validating deterministic Uniform-PAM baselines (%d physical points)...\n',nM*nSig*nSNR);
    for iM=1:nM
        M=MVec(iM); xPAM=define_constellation(M,P_avg,true,"mean");
        for iSig=1:nSig
            for iSNR=1:nSNR
                selectedDt=dtGrid(iM,iSig,iSNR); meta=dtMeta{iM,iSig,iSNR};
                cfg=build_fso_config(M,P_avg,snrVec(iSNR),sigVec(iSig), ...
                    'FastFadingMethod',o.FastFadingMethod,'ghN_h',o.ghN_h, ...
                    'loghDt',selectedDt,'loghSpanSigma',o.loghSpanSigma, ...
                    'loghYBlockSize',o.loghYBlockSize,'YGridMode',o.YGridMode, ...
                    'YGridScale',o.YGridScale,'xMaxBound',xMaxBounds(iM), ...
                    'saMaxIter',1,'saNStarts',1,'saUseParallel',false, ...
                    'minGap',o.minGap,'pinZero',o.pinZero,'seedInit',1,'logEvery',0,'historyEvery',1);
                AMI_functions.assert_constellation_feasible(xPAM,cfg,1e-10);
                t=tic; amiVal=cfg.AMI_Validator(xPAM(:)); valRuntime=toc(t); amiFast=cfg.AMI_Evaluator(xPAM(:));
                baseline{iM,iSig,iSNR}=struct('x',xPAM(:),'amiFast',amiFast, ...
                    'amiValidated',amiVal,'fastValidationGap',amiFast-amiVal, ...
                    'validationRuntime',valRuntime,'loghDtSelected',selectedDt, ...
                    'loghDtPolicy',meta.policy,'loghDtReason',meta.reason,'loghDtIsRefined',meta.isRefined);
            end
        end
    end
end


function [c,path]=run_revision_case(iM,iSig,iSNR,iR,MVec,sigVec,snrVec, ...
        P_avg,o,xMaxBounds,saProfile,baseline,runMeta,caseDir,dtGrid,dtMeta)
    M=MVec(iM); sigma_X_sq=sigVec(iSig); SNR_dB=snrVec(iSNR);
    selectedDt=dtGrid(iM,iSig,iSNR); meta=dtMeta{iM,iSig,iSNR};
    replicateBaseSeed=double(o.baseSeed)+(double(iR)-1)*10000019;
    seed=fso_result_utils.case_seed(replicateBaseSeed,M,SNR_dB,sigma_X_sq);

    cfg=build_fso_config(M,P_avg,SNR_dB,sigma_X_sq, ...
        'FastFadingMethod',o.FastFadingMethod,'ghN_h',o.ghN_h, ...
        'loghDt',selectedDt,'loghSpanSigma',o.loghSpanSigma, ...
        'loghYBlockSize',o.loghYBlockSize,'YGridMode',o.YGridMode,'YGridScale',o.YGridScale, ...
        'xMaxBound',xMaxBounds(iM),'saMaxIter',o.saMaxIter,'saNStarts',o.saNStarts, ...
        'saUseParallel',false,'minGap',o.minGap,'pinZero',o.pinZero, ...
        'seedInit',seed,'logEvery',0,'historyEvery',o.historyEvery);

    cfg.SA.T0=saProfile.T0; cfg.SA.Tf=saProfile.Tf; cfg.SA.baseStd0=saProfile.baseStd0;
    cfg.SA.itersPerTemp=saProfile.itersPerTemp; cfg.SA.nBlocks=ceil(cfg.SA.maxIter/cfg.SA.itersPerTemp);
    cfg.SA.coolingRate=exp(log(cfg.SA.Tf/cfg.SA.T0)/cfg.SA.nBlocks);

    xPAM=baseline.x(:); AMI_functions.assert_constellation_feasible(xPAM,cfg,1e-10);
    label=sprintf('[RevAMI|M%d|SNR%.1f|sig%.3f|rep%d|dt%.4g]',M,SNR_dB,sigma_X_sq,iR,selectedDt);
    tWall=tic; [saOut,starts]=sa_multistart(cfg,xPAM,label); wallRuntime=toc(tWall);

    win=starts(saOut.bestStartValidated);
    k95=convergence_iteration(win.history,win.initialMI,win.bestMIFast,0.95);
    k99=convergence_iteration(win.history,win.initialMI,win.bestMIFast,0.99);
    restartVals=[starts.bestMIValidated]; totalFastEvaluations=double(o.saNStarts)*double(o.saMaxIter);
    optimizerSaved=saOut; if isfield(optimizerSaved,'bestRun'),optimizerSaved=rmfield(optimizerSaved,'bestRun');end

    c=struct(); c.meta=runMeta;
    c.case=struct('M',M,'P_avg',P_avg,'SNR_dB',SNR_dB,'sigma_X_sq',sigma_X_sq,'replicate',iR,'seed',seed);
    c.integration=struct('loghDtSelected',selectedDt,'loghDtPolicy',meta.policy, ...
        'loghDtReason',meta.reason,'loghDtIsRefined',meta.isRefined,'loghDtRuleVersion',meta.ruleVersion);
    c.sa=struct('T0',cfg.SA.T0,'Tf',cfg.SA.Tf,'baseStd0',cfg.SA.baseStd0, ...
        'itersPerTemp',cfg.SA.itersPerTemp,'coolingRate',cfg.SA.coolingRate, ...
        'maxIter',cfg.SA.maxIter,'nStarts',cfg.SA.nStarts);
    c.config=fso_result_utils.config_snapshot(cfg); c.baseline=baseline; c.starts=starts; c.optimizer=optimizerSaved;
    c.winner=struct('startIndex',saOut.bestStartValidated,'x',saOut.bestXValidated(:), ...
        'amiFast',saOut.validatedWinnerFastMI,'amiValidated',saOut.bestMIValidated, ...
        'gainBits',saOut.bestMIValidated-baseline.amiValidated,'acceptanceRate',win.acceptanceRate, ...
        'convergenceIter95',k95,'convergenceIter99',k99,'history',win.history);
    c.repeatability=struct('restartValidatedMean',mean(restartVals,'omitnan'), ...
        'restartValidatedStd',std(restartVals,0,'omitnan'),'restartValidatedMin',min(restartVals,[],'omitnan'), ...
        'restartValidatedMax',max(restartVals,[],'omitnan'),'restartValidatedRange',max(restartVals)-min(restartVals));
    c.validation=struct('selectionChanged',saOut.selectionChanged,'fastWinnerIndex',saOut.bestStartFast, ...
        'validatedWinnerIndex',saOut.bestStartValidated,'validatedAdvantageOverFastWinner',saOut.validatedAdvantageOverFastWinner, ...
        'meanAbsFastValidationGap',saOut.meanAbsFastValidationGap,'maxAbsFastValidationGap',saOut.maxAbsFastValidationGap);
    c.geometry=struct('actualMinGap',safe_min_diff(saOut.bestXValidated), ...
        'xMax',max(saOut.bestXValidated),'meanPower',mean(saOut.bestXValidated));
    c.runtime=struct('wallSeconds',wallRuntime,'saSecondsSummed',saOut.totalSARuntime, ...
        'validationSeconds',saOut.validationTotalRuntime,'fastEvaluations',totalFastEvaluations, ...
        'fastEvaluationsPerSecond',totalFastEvaluations/max(saOut.totalSARuntime,eps));

    fname=sprintf('case_M%d_SNR%s_sig%s_rep%02d.mat',M,num_token(SNR_dB,2),num_token(sigma_X_sq,4),iR);
    path=fullfile(caseDir,fname); save(path,'c','-v7.3');
end


function S=build_summary(caseGrid,baseline,MVec,sigVec,snrVec,nRep)
    nM=numel(MVec); nSig=numel(sigVec); nSNR=numel(snrVec); sz4=[nM nSig nSNR nRep]; sz3=[nM nSig nSNR];
    S.raw=struct();
    fields={'amiGSValidated','amiGSFast','gainBits','convergenceIter95','convergenceIter99', ...
        'acceptanceRate','restartValidatedStd','meanAbsFastValidationGap','maxAbsFastValidationGap', ...
        'selectionChanged','wallSeconds','saSecondsSummed','validationSeconds','fastEvaluationsPerSecond', ...
        'actualMinGap','xMax','meanPower','loghDtSelected','loghDtIsRefined'};
    for k=1:numel(fields),S.raw.(fields{k})=nan(sz4);end
    S.xWinner=cell(sz4); S.amiPAMValidated=nan(sz3); S.amiPAMFast=nan(sz3); S.pamFastValidationGap=nan(sz3);
    S.loghDtSelected=nan(sz3); S.loghDtIsRefined=false(sz3);

    for iM=1:nM
        for iSig=1:nSig
            for iSNR=1:nSNR
                b=baseline{iM,iSig,iSNR};
                S.amiPAMValidated(iM,iSig,iSNR)=b.amiValidated; S.amiPAMFast(iM,iSig,iSNR)=b.amiFast;
                S.pamFastValidationGap(iM,iSig,iSNR)=b.fastValidationGap;
                S.loghDtSelected(iM,iSig,iSNR)=b.loghDtSelected; S.loghDtIsRefined(iM,iSig,iSNR)=b.loghDtIsRefined;
                for iR=1:nRep
                    c=caseGrid{iM,iSig,iSNR,iR};
                    S.raw.amiGSValidated(iM,iSig,iSNR,iR)=c.winner.amiValidated;
                    S.raw.amiGSFast(iM,iSig,iSNR,iR)=c.winner.amiFast;
                    S.raw.gainBits(iM,iSig,iSNR,iR)=c.winner.gainBits;
                    S.raw.convergenceIter95(iM,iSig,iSNR,iR)=c.winner.convergenceIter95;
                    S.raw.convergenceIter99(iM,iSig,iSNR,iR)=c.winner.convergenceIter99;
                    S.raw.acceptanceRate(iM,iSig,iSNR,iR)=c.winner.acceptanceRate;
                    S.raw.restartValidatedStd(iM,iSig,iSNR,iR)=c.repeatability.restartValidatedStd;
                    S.raw.meanAbsFastValidationGap(iM,iSig,iSNR,iR)=c.validation.meanAbsFastValidationGap;
                    S.raw.maxAbsFastValidationGap(iM,iSig,iSNR,iR)=c.validation.maxAbsFastValidationGap;
                    S.raw.selectionChanged(iM,iSig,iSNR,iR)=double(c.validation.selectionChanged);
                    S.raw.wallSeconds(iM,iSig,iSNR,iR)=c.runtime.wallSeconds;
                    S.raw.saSecondsSummed(iM,iSig,iSNR,iR)=c.runtime.saSecondsSummed;
                    S.raw.validationSeconds(iM,iSig,iSNR,iR)=c.runtime.validationSeconds;
                    S.raw.fastEvaluationsPerSecond(iM,iSig,iSNR,iR)=c.runtime.fastEvaluationsPerSecond;
                    S.raw.actualMinGap(iM,iSig,iSNR,iR)=c.geometry.actualMinGap;
                    S.raw.xMax(iM,iSig,iSNR,iR)=c.geometry.xMax; S.raw.meanPower(iM,iSig,iSNR,iR)=c.geometry.meanPower;
                    S.raw.loghDtSelected(iM,iSig,iSNR,iR)=c.integration.loghDtSelected;
                    S.raw.loghDtIsRefined(iM,iSig,iSNR,iR)=double(c.integration.loghDtIsRefined);
                    S.xWinner{iM,iSig,iSNR,iR}=c.winner.x;
                end
            end
        end
    end

    S.mean=struct(); S.stdAcrossReplicates=struct(); S.minAcrossReplicates=struct(); S.maxAcrossReplicates=struct();
    rawNames=fieldnames(S.raw);
    for k=1:numel(rawNames)
        f=rawNames{k}; a=S.raw.(f); S.mean.(f)=mean(a,4,'omitnan'); S.stdAcrossReplicates.(f)=std(a,0,4,'omitnan');
        S.minAcrossReplicates.(f)=min(a,[],4,'omitnan'); S.maxAcrossReplicates.(f)=max(a,[],4,'omitnan');
    end
    S.selectionChangedFraction=S.mean.selectionChanged;
    S.bestReplicateIndex=nan(sz3); S.bestX=cell(sz3); S.bestAMIValidated=nan(sz3); S.bestGainBits=nan(sz3);
    for iM=1:nM
        for iSig=1:nSig
            for iSNR=1:nSNR
                vals=reshape(S.raw.amiGSValidated(iM,iSig,iSNR,:),1,[]); [bestVal,bestR]=max(vals,[],'omitnan');
                S.bestReplicateIndex(iM,iSig,iSNR)=bestR; S.bestAMIValidated(iM,iSig,iSNR)=bestVal;
                S.bestGainBits(iM,iSig,iSNR)=S.raw.gainBits(iM,iSig,iSNR,bestR); S.bestX{iM,iSig,iSNR}=S.xWinner{iM,iSig,iSNR,bestR};
            end
        end
    end
end


function files=make_summary_plots(bundle,figDir,saveFigures)
    MVec=bundle.axes.M; sigVec=bundle.axes.sigma_X_sq; snrVec=bundle.axes.SNR_dB; S=bundle.summary; files={};
    for iSig=1:numel(sigVec)
        f=figure('Name',sprintf('Revised AMI sigma=%.3f',sigVec(iSig)),'Color','w');
        nCols=min(2,numel(MVec)); nRows=ceil(numel(MVec)/nCols); tl=tiledlayout(f,nRows,nCols,'TileSpacing','compact','Padding','compact');
        title(tl,sprintf('Validated AMI vs SNR, \\sigma_X^2=%.3f, d_{min}=%.3f',sigVec(iSig),bundle.method.minGap));
        for iM=1:numel(MVec)
            ax=nexttile(tl);hold(ax,'on');grid(ax,'on');box(ax,'on');
            yGS=reshape(S.mean.amiGSValidated(iM,iSig,:),1,[]); eGS=reshape(S.stdAcrossReplicates.amiGSValidated(iM,iSig,:),1,[]);
            yPAM=reshape(S.amiPAMValidated(iM,iSig,:),1,[]);
            errorbar(ax,snrVec,yGS,eGS,'-o','LineWidth',1.5,'DisplayName','GS mean +/- std');
            plot(ax,snrVec,yPAM,'--s','LineWidth',1.3,'DisplayName','Uniform PAM'); yline(ax,log2(MVec(iM)),':','HandleVisibility','off');
            xlabel(ax,'SNR [dB]');ylabel(ax,'AMI [bits/symbol]');title(ax,sprintf('M=%d',MVec(iM)));legend(ax,'Location','best');
        end
        if saveFigures,files=[files save_figure_pair(f,figDir,sprintf('AMI_vs_SNR_sig%s',num_token(sigVec(iSig),4)))];end %#ok<AGROW>
    end
    for iSig=1:numel(sigVec)
        f=figure('Name',sprintf('Revised shaping gain sigma=%.3f',sigVec(iSig)),'Color','w');ax=axes(f);hold(ax,'on');grid(ax,'on');box(ax,'on');
        for iM=1:numel(MVec)
            y=reshape(S.mean.gainBits(iM,iSig,:),1,[]);e=reshape(S.stdAcrossReplicates.gainBits(iM,iSig,:),1,[]);
            errorbar(ax,snrVec,y,e,'-o','LineWidth',1.5,'DisplayName',sprintf('M=%d',MVec(iM)));
        end
        xlabel(ax,'SNR [dB]');ylabel(ax,'Validated GS - Uniform PAM [bits/symbol]');
        title(ax,sprintf('Geometric-shaping gain, \\sigma_X^2=%.3f, d_{min}=%.3f',sigVec(iSig),bundle.method.minGap));legend(ax,'Location','best');
        if saveFigures,files=[files save_figure_pair(f,figDir,sprintf('shaping_gain_sig%s',num_token(sigVec(iSig),4)))];end %#ok<AGROW>
    end
    for iSig=1:numel(sigVec)
        f=figure('Name',sprintf('Revision diagnostics sigma=%.3f',sigVec(iSig)),'Color','w');
        tl=tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');title(tl,sprintf('Optimization diagnostics, \\sigma_X^2=%.3f',sigVec(iSig)));
        ax1=nexttile(tl,1);hold(ax1,'on');grid(ax1,'on');
        for iM=1:numel(MVec)
            y=reshape(S.mean.convergenceIter99(iM,iSig,:),1,[]);e=reshape(S.stdAcrossReplicates.convergenceIter99(iM,iSig,:),1,[]);
            errorbar(ax1,snrVec,y,e,'-o','LineWidth',1.3,'DisplayName',sprintf('M=%d',MVec(iM)));
        end
        yline(ax1,bundle.sa.maxIter,'--','Iteration budget');xlabel(ax1,'SNR [dB]');ylabel(ax1,'k_{99} [iterations]');title(ax1,'Convergence');legend(ax1,'Location','best');
        ax2=nexttile(tl,2);hold(ax2,'on');grid(ax2,'on');
        for iM=1:numel(MVec)
            y=reshape(S.stdAcrossReplicates.amiGSValidated(iM,iSig,:),1,[]);plot(ax2,snrVec,y,'-o','LineWidth',1.3,'DisplayName',sprintf('M=%d',MVec(iM)));
        end
        xlabel(ax2,'SNR [dB]');ylabel(ax2,'Std validated AMI [bits/symbol]');title(ax2,'Independent-run repeatability');legend(ax2,'Location','best');
        if saveFigures,files=[files save_figure_pair(f,figDir,sprintf('diagnostics_sig%s',num_token(sigVec(iSig),4)))];end %#ok<AGROW>
    end
end


function print_summary(bundle)
    S=bundle.summary; fprintf('\nREVISED AMI-vs-SNR SUMMARY\n');
    if strcmp(bundle.method.fastFadingMethod,'logh')
        if strcmp(bundle.method.loghDtPolicy,'case-adaptive')
            fprintf(['d_min=%.5g | fast=log-h(dt policy I9: base %.5g, refined %.5g; refined cases=%d) ' ...
                '| span=+/-%.3g sigma_t | y-grid=%s(scale=%.3g) | SA=%d x %d | replicates=%d\n'], ...
                bundle.method.minGap,bundle.method.loghDtBase,bundle.method.loghDtRefined, ...
                bundle.method.nRefinedPhysicalCases,bundle.method.loghSpanSigma, ...
                bundle.method.yGridMode,bundle.method.yGridScale,bundle.sa.nStarts,bundle.sa.maxIter,bundle.sa.nReplicates);
        else
            fprintf(['d_min=%.5g | fast=log-h(fixed dt=%.5g, span=+/-%.3g sigma_t) ' ...
                '| y-grid=%s(scale=%.3g) | SA=%d x %d | replicates=%d\n'], ...
                bundle.method.minGap,bundle.method.loghDt,bundle.method.loghSpanSigma, ...
                bundle.method.yGridMode,bundle.method.yGridScale,bundle.sa.nStarts,bundle.sa.maxIter,bundle.sa.nReplicates);
        end
    else
        fprintf('d_min=%.5g | fast=GH%d | y-grid=fixed-bound | SA=%d x %d | replicates=%d\n', ...
            bundle.method.minGap,bundle.method.ghN_h,bundle.sa.nStarts,bundle.sa.maxIter,bundle.sa.nReplicates);
    end
    fprintf('SA profile: T0=%.5g Tf=%.5g step0=%.5g block=%d\n',bundle.sa.T0,bundle.sa.Tf,bundle.sa.baseStd0,bundle.sa.itersPerTemp);
    for iSig=1:numel(bundle.axes.sigma_X_sq)
        fprintf('\n--- sigma_X^2 = %.4f ---\n',bundle.axes.sigma_X_sq(iSig));
        for iM=1:numel(bundle.axes.M)
            ami=reshape(S.mean.amiGSValidated(iM,iSig,:),1,[]); sd=reshape(S.stdAcrossReplicates.amiGSValidated(iM,iSig,:),1,[]);
            gain=reshape(S.mean.gainBits(iM,iSig,:),1,[]); k99=reshape(S.mean.convergenceIter99(iM,iSig,:),1,[]);
            fv=reshape(S.mean.maxAbsFastValidationGap(iM,iSig,:),1,[]); dt=reshape(S.loghDtSelected(iM,iSig,:),1,[]);
            fprintf('M=%d\n',bundle.axes.M(iM)); fprintf('  SNR       : %s\n',mat2str(bundle.axes.SNR_dB));
            fprintf('  logh dt   : %s\n',mat2str(dt,5)); fprintf('  GS mean   : %s\n',mat2str(ami,6));
            fprintf('  GS std    : %s\n',mat2str(sd,4)); fprintf('  gain mean : %s\n',mat2str(gain,6));
            fprintf('  k99 mean  : %s\n',mat2str(k99,6)); fprintf('  F-V gap   : %s\n',mat2str(fv,4));
        end
    end
    fprintf('\nGlobal diagnostics:\n');
    fprintf('  max between-replicate AMI std: %.3e bits/symbol\n',max(S.stdAcrossReplicates.amiGSValidated(:),[],'omitnan'));
    fprintf('  max fast-validation restart gap: %.3e bits/symbol\n',max(S.raw.maxAbsFastValidationGap(:),[],'omitnan'));
    fprintf('  selection-changed fraction: %.3f\n',mean(S.raw.selectionChanged(:),'omitnan'));
    fprintf('  max k99 observed: %.0f / %.0f iterations\n',max(S.raw.convergenceIter99(:),[],'omitnan'),bundle.sa.maxIter);
end

function k=convergence_iteration(history,initialMI,bestMI,fraction)
    improvement=bestMI-initialMI;if ~(isfinite(improvement)&&improvement>0),k=0;return;end
    target=initialMI+fraction*improvement;idx=find(history.bestMI>=target,1,'first');if isempty(idx),k=history.iter(end);else,k=history.iter(idx);end
end
function validate_spacing_feasibility(MVec,P_avg,minGap)
    required=(MVec-1)*minGap/2;bad=find(required>P_avg+100*eps(max(1,P_avg)),1,'first');
    if ~isempty(bad),error('sim_revision_AMI_vs_SNR:InfeasibleGap','M=%d with d_min=%.6g requires mean power %.6g > P_avg=%.6g.',MVec(bad),minGap,required(bad),P_avg);end
end
function d=safe_min_diff(x),x=sort(real(x(:)),'ascend');if numel(x)<2,d=NaN;else,d=min(diff(x));end;end
function files=save_figure_pair(fig,figDir,baseName),if ~exist(figDir,'dir'),mkdir(figDir);end;pngPath=fullfile(figDir,[baseName '.png']);figPath=fullfile(figDir,[baseName '.fig']);exportgraphics(fig,pngPath,'Resolution',220);savefig(fig,figPath);files={pngPath,figPath};end
function token=sanitize_run_label(label),token=strtrim(char(string(label)));if isempty(token),token='run';return;end;token=regexprep(token,'[^A-Za-z0-9._-]+','_');token=regexprep(token,'_+','_');token=regexprep(token,'^[._-]+|[._-]+$','');if isempty(token),token='run';end;end
function [outDir,runId]=make_unique_run_directory(parentDir,requestedRunId),if ~exist(parentDir,'dir'),mkdir(parentDir);end;runId=requestedRunId;outDir=fullfile(parentDir,runId);suffix=1;while exist(outDir,'dir'),suffix=suffix+1;runId=sprintf('%s_%02d',requestedRunId,suffix);outDir=fullfile(parentDir,runId);end;[ok,msg]=mkdir(outDir);if ~ok,error('sim_revision_AMI_vs_SNR:RunDirectoryCreateFailed','Could not create run directory %s: %s',outDir,msg);end;end
function token=num_token(v,nDec),token=sprintf(['%0.' num2str(nDec) 'f'],double(v));token=strrep(token,'-','m');token=strrep(token,'+','p');token=strrep(token,'.','p');end
function tf=positive_scalar(v),tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&isreal(v)&&v>0;end
function tf=positive_integer(v),tf=positive_scalar(v)&&mod(v,1)==0;end
