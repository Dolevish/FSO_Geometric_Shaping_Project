function R5 = fig7_m8_20db_sigma01(varargin)
%FIG7_M8_20DB_SIGMA01 Commit K5: the missing publication Fig.7 comparison.
%
% Exact I.13D Pattern Search scientific policy at M=8, SNR=20 dB,
% sigma_R^2=0.1. Reuses frozen I.11.1 PAM/SA and I.13B Zhou results.
% Only Pattern Search is evaluated: 2 replicates x 6 restarts,
% max 4001 objective evaluations/start, independent validation/start.
% No parameter tuning, no SA/Zhou re-optimization, no global-optimum claim.
%
% dry=fig7_m8_20db_sigma01('DryRun',true);
% R5=fig7_m8_20db_sigma01();
% R5=fig7_m8_20db_sigma01('ResumeDirectory','/actual/K5_...');

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'ProductionBundlePath','',@text_scalar);
    addParameter(p,'ZhouResultPath','',@text_scalar);
    addParameter(p,'OutputDirectory','',@text_scalar);
    addParameter(p,'ResumeDirectory','',@text_scalar);
    addParameter(p,'ParallelWorkers',8,@positive_integer);
    addParameter(p,'FreshPool',true,@logical_scalar);
    addParameter(p,'DeletePoolOnFinish',true,@logical_scalar);
    addParameter(p,'DryRun',false,@logical_scalar);
    parse(p,varargin{:}); o=p.Results;

    original=i13d_benchmark_config();
    P=i11p1_production_config();
    % Abort if the I.13D scientific policy has changed. The selected
    % condition is NEW; the optimizer/constraints/numerics are NOT.
    assert(isequal(original.MVec,[4 8 16 32]) && ...
        isequal(original.SNRVec,[20 30]) && ...
        isequal(original.TurbulenceVec,[0 0.3]) && ...
        isequal(original.Replicates,[1 2]), ...
        'K5: unexpected I.13D benchmark grid.');
    assert(original.P_avg==1 && original.minGap==0.01 && ...
        original.pinZero && P.P_avg==original.P_avg && ...
        P.minGap==original.minGap && P.pinZero==original.pinZero, ...
        'K5: frozen mean-power, zero-anchor or minimum-gap policy changed.');
    assert(isequal(original.RestartsByM,[6 6 8 8]) && ...
        isequal(original.MaxIterByM,[4000 4000 4000 8000]) && ...
        isequal(original.MaxEvalsPerStart,[4001 4001 4001 8001]), ...
        'K5: frozen I.13D restart/evaluation budgets changed.');
    assert(strcmp(P.FastFadingMethod,'logh') && ...
        strcmp(P.YGridMode,'candidate-adaptive') && ...
        P.YGridScale==1.1 && P.loghSpanSigma==8 && ...
        P.loghYBlockSize==512 && P.historyEvery==50 && ...
        strcmp(original.PatternSearch.Display,'off') && ...
        ~original.PatternSearch.UseCompletePoll && ...
        ~original.PatternSearch.UseParallel && ...
        original.PatternSearch.MeshTolerance==1e-6 && ...
        original.PatternSearch.StepTolerance==1e-8 && ...
        original.PatternSearch.FunctionTolerance==1e-8 && ...
        original.PatternSearch.MaxFinalRepairInfNorm==1e-5, ...
        'K5: frozen evaluator or Pattern Search policy changed.');

    F=original;
    F.version='K5-FIG7-I13D-v1';
    % Subselect ONE physical point only. All original optimizer options
    % and constraints are inherited, not redefined or retuned.
    F.MVec=8; F.SNRVec=20; F.TurbulenceVec=0.1;
    F.RestartsByM=original.RestartsByM(original.MVec==8);
    F.MaxIterByM=original.MaxIterByM(original.MVec==8);
    F.MaxEvalsPerStart=original.MaxEvalsPerStart(original.MVec==8);
    F.nPhysicalPoints=1; F.nOptimizationCases=2;
    F.nRestartTasks=2*F.RestartsByM;
    F.totalMaxFunctionEvaluations=F.nRestartTasks*F.MaxEvalsPerStart;

    prodPath=resolve_production(o.ProductionBundlePath);
    X=load(prodPath,'bundle');
    if ~isfield(X,'bundle'),error('K5:Production','Missing production bundle.');end
    B=X.bundle;
    if ~isfield(B,'productionFreeze') || ...
            ~strcmp(B.productionFreeze.freezeVersion,P.freezeVersion)
        error('K5:ProductionFreeze','Expected original frozen I.11.1 production.');
    end
    [casePlan,taskPlan]=i13d_build_task_plan(B.summary.caseTable,F);
    if height(casePlan)~=2 || height(taskPlan)~=12 || ...
            ~isequal(sort(casePlan.Replicate).',[1 2]) || ...
            any(casePlan.NStarts~=6) || any(casePlan.MaxIterSA~=4000) || ...
            any(casePlan.MaxEvalsPerStart~=4001) || ...
            any(abs(casePlan.SelectedDt-0.01)>1e-12)
        error('K5:Plan','Unexpected I.11.1 production plan for M8/20/0.1.');
    end

    zhouPath=resolve_zhou(o.ZhouResultPath,prodPath);
    Z=load(zhouPath,'result');
    if ~isfield(Z,'result') || ~strcmp(Z.result.version,'I13B-v1') || ...
            ~strcmp(canonical(Z.result.sourceProduction),canonical(prodPath))
        error('K5:ZhouSource', ...
            'The I.13B Zhou result must come from the SAME production run.');
    end
    Q=Z.result.comparison;
    q=Q.M==8 & abs(Q.SNRdB-20)<1e-12 & abs(Q.SigmaX2-0.1)<1e-12;
    if nnz(q)~=1,error('K5:ZhouRow','Expected exactly one Zhou row.');end
    pam=mean(casePlan.PAMValidated);
    sa=mean(casePlan.SAValidated);
    zhou=Q.ZhouValidated(q);
    if ~all(isfinite([pam sa zhou])) || ...
            abs(Q.PAMValidated(q)-pam)>1e-8 || ...
            abs(Q.GSMeanValidated(q)-sa)>1e-8
        error('K5:SourceMismatch','Independent I.13B vs I.11.1 scores disagree.');
    end

    stage=B.stages.standard;
    iM=find(stage.axes.M==8,1);
    iT=find(abs(stage.axes.sigma_X_sq-0.1)<1e-12,1);
    iS=find(abs(stage.axes.SNR_dB-20)<1e-12,1);
    if numel(iM)~=1 || numel(iT)~=1 || numel(iS)~=1
        error('K5:Axes','Frozen production channel axes are not unique.');
    end
    base=stage.baseline{iM,iT,iS};
    if abs(base.amiValidated-pam)>1e-8
        error('K5:PAMBaseline','The original PAM baseline differs.');
    end

    specs=cell(2,1);
    for c=1:2
        row=casePlan(c,:);
        % Identical arguments/order to build_case_specs in I.13D.
        cfg=build_fso_config(row.M,F.P_avg,row.SNRdB,row.SigmaR2, ...
            'FastFadingMethod',P.FastFadingMethod, ...
            'loghDt',row.SelectedDt, ...
            'loghSpanSigma',P.loghSpanSigma, ...
            'loghYBlockSize',P.loghYBlockSize, ...
            'YGridMode',P.YGridMode,'YGridScale',P.YGridScale, ...
            'minGap',F.minGap,'pinZero',F.pinZero, ...
            'saMaxIter',row.MaxIterSA,'saNStarts',row.NStarts, ...
            'seedInit',row.Seed,'logEvery',0,'historyEvery',P.historyEvery);
        AMI_functions.assert_constellation_feasible(base.x(:),cfg,1e-10);
        if cfg.M~=8 || cfg.P_avg~=1 || cfg.SA.minGap~=0.01 || ...
                ~cfg.SA.pinZero || cfg.SA.maxIter~=4000 || ...
                cfg.SA.nStarts~=6 || abs(cfg.loghDt-0.01)>1e-12 || ...
                ~strcmp(cfg.fastFadingMethod,'logh') || ...
                ~strcmp(cfg.yGridMode,'candidate-adaptive')
            error('K5:ScientificConfig','Unexpected fast evaluator/constraint config.');
        end
        specs{c}=struct('cfg',cfg,'xBaseline',base.x(:), ...
            'seed',uint32(row.Seed),'maxEvals',row.MaxEvalsPerStart);
    end

    signature=struct('version',F.version, ...
        'sourceProduction',canonical(prodPath),'sourceZhou',canonical(zhouPath), ...
        'scientificConfig',P,'benchmarkConfig',F, ...
        'caseIdentities',table2struct(casePlan(:, ...
        {'M','SNRdB','SigmaR2','Replicate','Seed','SelectedDt', ...
        'NStarts','MaxIterSA','MaxEvalsPerStart'})));
    sigText=char(jsonencode(signature));

    fprintf('\nCOMMIT K5 Fig.7 | M=8 SNR=20 dB sigma_R^2=0.1\n');
    fprintf('Original validated PAM=%.8f | Zhou=%.8f | SA mean=%.8f\n', ...
        pam,zhou,sa);
    fprintf('Frozen I.13D: 2 replicates x 6 starts x 4001 evaluations = 12 tasks\n');
    fprintf('P_avg=1 | pinZero=1 | d_min=0.01 | logh dt=0.01\n');
    fprintf('All source files and all scientific settings are read-only.\n');
    if o.DryRun
        R5=struct('dryRun',true,'config',F,'casePlan',casePlan, ...
            'taskPlan',taskPlan,'signature',signature, ...
            'PAM',pam,'Zhou',zhou,'SA',sa,'outputDirectory','');
        return;
    end
    if exist('patternsearch','file')~=2
        error('K5:Toolbox','Global Optimization Toolbox is required.');
    end

    outDir=prepare_directory(o);
    checkpointDir=fullfile(outDir,'checkpoints');
    if exist(checkpointDir,'dir')~=7,mkdir(checkpointDir);end
    manifestPath=fullfile(outDir,'K5_scientific_manifest.mat');
    if exist(manifestPath,'file')==2
        S=load(manifestPath,'sigText');
        if ~isfield(S,'sigText') || ~strcmp(S.sigText,sigText)
            error('K5:ResumeSignature','Saved run has a different scientific signature.');
        end
    else
        if strlength(string(o.ResumeDirectory))>0
            error('K5:ResumeManifest','Missing K5 scientific manifest.');
        end
        save(manifestPath,'sigText','signature');
    end
    fprintf('K5 run directory: %s\n',outDir);

    n=height(taskPlan);
    tasks=cell(n,1);
    for j=1:n
        [~,~,identity,path]=task_args(j,taskPlan,casePlan,checkpointDir);
        if exist(path,'file')~=2,continue;end
        S=load(path,'checkpoint');
        if ~isfield(S,'checkpoint') || ...
                ~isfield(S.checkpoint,'result') || ...
                ~strcmp(S.checkpoint.signatureText,sigText) || ...
                ~isequaln(S.checkpoint.identity,identity) || ...
                ~S.checkpoint.result.success
            error('K5:Checkpoint','Invalid checkpoint: %s',path);
        end
        tasks{j}=S.checkpoint.result;
    end
    pending=find(cellfun(@isempty,tasks));
    fprintf('Reused validated checkpoints %d/12 | pending %d\n',n-numel(pending),numel(pending));

    if ~isempty(pending)
        if o.ParallelWorkers==1
            for i=1:numel(pending)
                j=pending(i);
                [c,t,id,path]=task_args(j,taskPlan,casePlan,checkpointDir);
                s=specs{c};
                tasks{j}=i13d_patternsearch_restart_task(s.cfg, ...
                    s.xBaseline,s.seed,t.Restart,s.maxEvals, ...
                    F.PatternSearch,id,path,sigText);
                validate_task(tasks{j},casePlan(c,:));
                fprintf('K5 DONE %d/%d | %s\n',i,numel(pending),id.IdentityKey);
            end
        else
            pp=gcp('nocreate');
            if ~isempty(pp)&&o.FreshPool,delete(pp);pp=[];end
            if isempty(pp)
                pp=parpool('Processes',double(o.ParallelWorkers));
            elseif pp.NumWorkers~=o.ParallelWorkers
                delete(pp);pp=parpool('Processes',double(o.ParallelWorkers));
            end
            if o.DeletePoolOnFinish
                cleanup=onCleanup(@() delete_pool(pp)); %#ok<NASGU>
            end
            active=parallel.FevalFuture.empty(0,1);
            activeRows=zeros(0,1);
            next=1;done=0;
            try
                while numel(active)<pp.NumWorkers && next<=numel(pending)
                    j=pending(next);
                    active(end+1,1)=submit_one(pp,j,taskPlan, ...
                        casePlan,specs,F,checkpointDir,sigText); %#ok<AGROW>
                    activeRows(end+1)=j; %#ok<AGROW>
                    next=next+1;
                end
                while ~isempty(active)
                    [k,res]=fetchNext(active);
                    j=activeRows(k);
                    active(k)=[];activeRows(k)=[];
                    c=taskPlan.TaskCaseIndex(j);
                    validate_task(res,casePlan(c,:));
                    tasks{j}=res;done=done+1;
                    fprintf('K5 DONE %d/%d | %s | validated=%.8f\n', ...
                        done,numel(pending),taskPlan.IdentityKey(j),res.bestMIValidated);
                    if next<=numel(pending)
                        j2=pending(next);
                        active(end+1,1)=submit_one(pp,j2,taskPlan, ...
                            casePlan,specs,F,checkpointDir,sigText); %#ok<AGROW>
                        activeRows(end+1)=j2; %#ok<AGROW>
                        next=next+1;
                    end
                end
            catch ME
                if ~isempty(active),try,cancel(active);catch,end,end
                rethrow(ME);
            end
        end
    end

    Replicate=zeros(n,1);Restart=zeros(n,1);PSFast=zeros(n,1);
    PSValidated=zeros(n,1);FunctionEvaluations=zeros(n,1);
    ValidationGap=zeros(n,1);TaskWallSeconds=zeros(n,1);
    for j=1:n
        r=tasks{j};
        c=taskPlan.TaskCaseIndex(j);
        validate_task(r,casePlan(c,:));
        Replicate(j)=r.Replicate;Restart(j)=r.Restart;
        PSFast(j)=r.bestMIFast;PSValidated(j)=r.bestMIValidated;
        FunctionEvaluations(j)=r.functionEvaluations;
        ValidationGap(j)=r.validationGap;TaskWallSeconds(j)=r.taskWallSeconds;
    end
    restartSummary=sortrows(table(Replicate,Restart,PSFast,PSValidated, ...
        FunctionEvaluations,ValidationGap,TaskWallSeconds), ...
        {'Replicate','Restart'});

    PSWinner=zeros(2,1); WinnerRestart=zeros(2,1);
    for rep=1:2
        A=restartSummary(restartSummary.Replicate==rep,:);
        if height(A)~=6 || ~isequal(sort(A.Restart).',1:6)
            error('K5:Replicate','Missing or duplicate restart.');
        end
        [PSWinner(rep),ii]=max(A.PSValidated);
        WinnerRestart(rep)=A.Restart(ii);
    end
    ps=mean(PSWinner);
    methods=["Uniform PAM";"Zhou GS";"Pattern Search";"SA"];
    vals=[pam;zhou;ps;sa];
    stds=[0;0;std(PSWinner,0);std(casePlan.SAValidated,0)];
    fourMethodTable=table(methods,vals,stds, ...
        'VariableNames',{'Method','ValidatedAMI','StdAcrossReplicates'});
    writetable(restartSummary,fullfile(outDir,'Fig7_patternsearch_restart_results.csv'));
    writetable(fourMethodTable,fullfile(outDir,'Fig7_four_method_values.csv'));
    figPaths=plot_fig7_four_methods(fourMethodTable,0.1,outDir);

    R5=struct('version',F.version,'config',F, ...
        'sourceProduction',prodPath,'sourceZhou',zhouPath, ...
        'casePlan',casePlan,'taskPlan',taskPlan, ...
        'restartSummary',restartSummary, ...
        'replicatePatternSearchAMI',PSWinner, ...
        'replicateWinnerRestart',WinnerRestart, ...
        'fourMethodTable',fourMethodTable,'figurePaths',figPaths, ...
        'outputDirectory',outDir,'globalOptimumClaim',false);
    save(fullfile(outDir,'Fig7_four_method_results.mat'),'R5','-v7.3');
    fprintf('\nCOMMIT K5 COMPLETE: independently validated four-method Fig.7\n');
    disp(fourMethodTable);
    fprintf('Fig.7 exported to %s\n',outDir);
end


function fut=submit_one(pp,j,taskPlan,casePlan,specs,F,checkpointDir,sigText)
    [c,t,id,path]=task_args(j,taskPlan,casePlan,checkpointDir);
    s=specs{c};
    fut=parfeval(pp,@i13d_patternsearch_restart_task,1, ...
        s.cfg,s.xBaseline,s.seed,t.Restart,s.maxEvals, ...
        F.PatternSearch,id,path,sigText);
end

function [c,t,id,path]=task_args(j,taskPlan,casePlan,checkpointDir)
    t=taskPlan(j,:);
    c=t.TaskCaseIndex;
    row=casePlan(c,:);
    id=struct('IdentityKey',char(t.IdentityKey), ...
        'M',row.M,'SNRdB',row.SNRdB,'SigmaR2',row.SigmaR2, ...
        'Replicate',row.Replicate,'Restart',t.Restart,'Seed',double(row.Seed));
    path=fullfile(checkpointDir,[char(t.IdentityKey) '.mat']);
end

function validate_task(r,C)
    if isempty(r)||~isstruct(r)||~isfield(r,'success')||~r.success
        if isstruct(r)&&isfield(r,'errorReport')
            error('K5:Worker','Pattern Search failed: %s',char(r.errorReport));
        end
        error('K5:Worker','Pattern Search restart failed.');
    end
    if ~isfinite(r.bestMIValidated)||~isfinite(r.bestMIFast)|| ...
            ~isfinite(r.functionEvaluations)||r.functionEvaluations>4001 || ...
            r.functionEvaluations<1 || ...
            r.M~=C.M || r.SNRdB~=C.SNRdB || ...
            abs(r.SigmaR2-C.SigmaR2)>1e-12 || ...
            r.Replicate~=C.Replicate || ...
            norm(r.xBest(:),1)<=0
        error('K5:WorkerValidity','Validation, identity, or PS budget mismatch.');
    end
end

function path=resolve_production(given)
    if strlength(string(given))>0
        path=char(string(given));
    else
        root=fullfile(fso_result_utils.default_results_root(),'i11p1_production');
        D=dir(fullfile(root,'*','i11p1ProductionBundle.mat'));
        if numel(D)~=1
            error('K5:ProductionDiscovery', ...
                'Expected one production bundle, found %d. Supply ProductionBundlePath explicitly.',numel(D));
        end
        path=fullfile(D(1).folder,D(1).name);
    end
    if exist(path,'file')~=2
        error('K5:ProductionMissing','Missing production MAT: %s',path);
    end
end

function path=resolve_zhou(given,prodPath)
    if strlength(string(given))>0
        path=char(string(given));
    else
        path=fullfile(fileparts(prodPath),'i13_reviewer12', ...
            'zhou_baseline','i13b_zhou_baseline.mat');
    end
    if exist(path,'file')~=2
        error('K5:ZhouMissing','Missing Zhou baseline MAT: %s',path);
    end
end

function outDir=prepare_directory(o)
    root=char(java.io.File(fso_result_utils.revision2_results_root()).getCanonicalPath());
    if strlength(string(o.ResumeDirectory))>0
        outDir=char(string(o.ResumeDirectory));
        if exist(outDir,'dir')~=7
            error('K5:ResumeMissing','Missing resume directory.');
        end
    elseif strlength(string(o.OutputDirectory))>0
        if exist(char(string(o.OutputDirectory)),'dir')==7
            error('K5:OutputExists','New K5 output directory must not already exist; use ResumeDirectory.');
        end
        outDir=char(string(o.OutputDirectory));
    else
        parent=fullfile(root,'commit_k5_fig7');
        if exist(parent,'dir')~=7,mkdir(parent);end
        token=char(datetime('now','TimeZone','UTC', ...
            'Format','yyyyMMdd''T''HHmmssSSS''Z'''));
        outDir=fullfile(parent,['K5_' token]);
    end
    target=canonical(outDir);
    if ~(strcmp(target,root)||startsWith(target,[root filesep]))
        error('K5:OutputRoot','Outputs must remain under results/ieee_revision2.');
    end
    if exist(outDir,'dir')~=7
        [ok,msg]=mkdir(outDir);
        if ~ok,error('K5:mkdir','%s',msg);end
    end
end

function s=canonical(path)
    s=char(java.io.File(char(string(path))).getCanonicalPath());
end
function delete_pool(pp)
    if ~isempty(pp),try,delete(pp);catch,end,end
end
function tf=text_scalar(v),tf=ischar(v)||(isstring(v)&&isscalar(v));end
function tf=positive_integer(v)
    tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&v==floor(v)&&v>=1;
end
function tf=logical_scalar(v),tf=islogical(v)&&isscalar(v);end
