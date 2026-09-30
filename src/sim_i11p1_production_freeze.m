function bundle = sim_i11p1_production_freeze(varargin)
%SIM_I11P1_PRODUCTION_FREEZE  Canonical final-production entry point.
%
% Frozen scientific policy:
%   M=4,8,16 -> 4000 iterations/restart
%   M=32     -> 8000 iterations/restart
%   restarts -> [6 6 8 8] for M=[4 8 16 32]
%
% Execution is intentionally staged: M=32 runs first with the full flat
% restart-task pool, followed by M=4,8,16. Each stage uses the already
% validated I.11 checkpointed flat scheduler and therefore preserves
% deterministic Threefry restart streams and per-restart validation.
%
% New run:
%   b = sim_i11p1_production_freeze();
%
% Resume the whole production run:
%   b = sim_i11p1_production_freeze('ResumeDirectory','/path/to/I11p1_run');

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'UseParallel',true,@logical_scalar);
    addParameter(p,'ParallelWorkers',8,@(v) isempty(v)||positive_integer(v));
    addParameter(p,'LongestFirst',true,@logical_scalar);
    addParameter(p,'WriteCSV',true,@logical_scalar);
    addParameter(p,'RunLabel','I11p1_production',@text_scalar);
    addParameter(p,'ResultsRoot',fso_result_utils.default_results_root(),@text_scalar);
    addParameter(p,'ResumeDirectory','',@text_scalar);
    addParameter(p,'DryRun',false,@logical_scalar);
    parse(p,varargin{:}); o=p.Results;

    F=i11p1_production_config();
    stages=i11p1_stage_plan(F);

    fprintf('\n============================================================\n');
    fprintf('I.11.1 PRODUCTION FREEZE (%s)\n',F.freezeVersion);
    fprintf('Grid: M=%s | SNR=%s dB | sigma_X^2=%s | replicates=%d\n', ...
        mat2str(F.MVec),mat2str(F.SNRVec),mat2str(F.TurbulenceVec),F.nReplicates);
    fprintf('Restarts by M   : %s\n',mat2str(F.RestartsByM));
    fprintf('Iterations by M : %s\n',mat2str(F.MaxIterByM));
    fprintf('Optimization cases=%d | restart tasks=%d | SA iterations=%d\n', ...
        F.nOptimizationCases,F.totalRestartTasks,F.totalSAIterations);
    fprintf('Stage order: M32 first (8000), then M4/8/16 (4000)\n');
    fprintf('No global-optimum claim is implied by this freeze.\n');
    fprintf('============================================================\n');

    if o.DryRun
        bundle=struct('productionFreeze',F,'stagePlan',stages,'dryRun',true);
        return;
    end

    signatureText=char(jsonencode(F));
    resumeRequested=strlength(string(o.ResumeDirectory))>0;
    if resumeRequested
        outDir=char(string(o.ResumeDirectory));
        manifestPath=fullfile(outDir,'production_manifest.mat');
        if exist(manifestPath,'file')~=2
            error('sim_i11p1_production_freeze:ResumeManifestMissing', ...
                'Resume directory has no production_manifest.mat: %s',outDir);
        end
        S=load(manifestPath,'manifest'); manifest=S.manifest;
        if ~isfield(manifest,'signatureText') || ~strcmp(manifest.signatureText,signatureText)
            error('sim_i11p1_production_freeze:ResumeConfigMismatch', ...
                'Resume directory does not match the frozen I.11.1 scientific profile.');
        end
    else
        parent=fullfile(char(o.ResultsRoot),'i11p1_production');
        if exist(parent,'dir')~=7,mkdir(parent);end
        token=char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmssSSS''Z'''));
        runId=sprintf('%s_%s',sanitize_label(o.RunLabel),token);
        outDir=fullfile(parent,runId); k=1;
        while exist(outDir,'dir')==7,k=k+1;outDir=fullfile(parent,sprintf('%s_%02d',runId,k));end
        [ok,msg]=mkdir(outDir); if ~ok,error('sim_i11p1_production_freeze:mkdir','%s',msg);end
        manifestPath=fullfile(outDir,'production_manifest.mat');
        manifest=struct();
        manifest.freezeVersion=F.freezeVersion;
        manifest.signatureText=signatureText;
        manifest.stagePaths=struct('m32','','standard','');
        manifest.stageCompleted=struct('m32',false,'standard',false);
        manifest.createdUTC=char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss.SSS''Z'''));
        save_manifest_atomic(manifestPath,manifest);
    end

    fprintf('Production root: %s\n',outDir);

    if o.UseParallel && ~isempty(o.ParallelWorkers)
        try
            pp=gcp('nocreate');
            if ~isempty(pp) && pp.NumWorkers~=double(o.ParallelWorkers),delete(pp);end
        catch
        end
    end

    stageBundles=struct();
    for k=1:numel(stages)
        st=stages(k); field=st.name;
        stageRoot=fullfile(outDir,[field '_stage']);
        if exist(stageRoot,'dir')~=7,mkdir(stageRoot);end

        childPath=manifest.stagePaths.(field);
        if isempty(childPath)
            childPath=discover_existing_child(stageRoot);
            if ~isempty(childPath)
                manifest.stagePaths.(field)=childPath;
                save_manifest_atomic(manifestPath,manifest);
            end
        end

        bundlePath='';
        if ~isempty(childPath)
            bundlePath=fullfile(childPath,'revisedAMIvsSNR_checkpointed_flat.mat');
        end
        if manifest.stageCompleted.(field) && exist(bundlePath,'file')==2
            L=load(bundlePath,'bundle'); b=L.bundle;
            fprintf('Stage %s already complete; reusing saved stage bundle.\n',field);
        else
            args=common_stage_args(F,st,o,stageRoot,field);
            if ~isempty(childPath)
                b=sim_revision_AMI_vs_SNR_checkpointed_flat(args{:},'ResumeDirectory',childPath);
            else
                b=sim_revision_AMI_vs_SNR_checkpointed_flat(args{:});
            end
            manifest.stagePaths.(field)=b.run.outputDirectory;
            manifest.stageCompleted.(field)=true;
            save_manifest_atomic(manifestPath,manifest);
        end
        stageBundles.(field)=b;
    end

    T=[stageBundles.standard.summary.caseTable; stageBundles.m32.summary.caseTable];
    T=sortrows(T,{'M','SigmaX2','SNRdB','Replicate'});

    bundle=struct();
    bundle.productionFreeze=F;
    bundle.stagePlan=stages;
    bundle.run=struct('outputDirectory',outDir,'resumeRequested',resumeRequested, ...
        'manifestPath',manifestPath);
    bundle.stages=stageBundles;
    bundle.summary=struct('caseTable',T);
    bundle.checkpoint=struct( ...
        'nReused',stageBundles.m32.checkpoint.nReused+stageBundles.standard.checkpoint.nReused, ...
        'nExecuted',stageBundles.m32.checkpoint.nExecuted+stageBundles.standard.checkpoint.nExecuted);

    if o.WriteCSV
        writetable(T,fullfile(outDir,'i11p1_combined_case_summary.csv'));
    end
    save(fullfile(outDir,'i11p1ProductionBundle.mat'),'bundle','-v7.3');
    fprintf('Saved canonical I.11.1 production bundle: %s\n', ...
        fullfile(outDir,'i11p1ProductionBundle.mat'));
end

function args=common_stage_args(F,st,o,stageRoot,field)
    args={'MVec',st.MVec,'P_avg',F.P_avg,'SNRVec',F.SNRVec, ...
        'TurbulenceVec',F.TurbulenceVec,'minGap',F.minGap,'pinZero',F.pinZero, ...
        'FastFadingMethod',F.FastFadingMethod,'LoghDtPolicy',F.LoghDtPolicy, ...
        'loghDt',F.loghDt,'loghDtRefined',F.loghDtRefined, ...
        'loghSpanSigma',F.loghSpanSigma,'loghYBlockSize',F.loghYBlockSize, ...
        'YGridMode',F.YGridMode,'YGridScale',F.YGridScale, ...
        'saMaxIter',st.saMaxIter,'RestartPolicy',F.RestartPolicy, ...
        'nReplicates',F.nReplicates,'T0',F.T0,'Tf',F.Tf, ...
        'BaseStd0',F.BaseStd0,'ItersPerTemp',F.ItersPerTemp, ...
        'baseSeed',F.baseSeed,'historyEvery',F.historyEvery, ...
        'UseParallel',o.UseParallel,'ParallelWorkers',o.ParallelWorkers, ...
        'LongestFirst',o.LongestFirst,'Checkpoint',true,'WriteCSV',o.WriteCSV, ...
        'RunLabel',sprintf('%s_%s',char(string(o.RunLabel)),field),'ResultsRoot',stageRoot};
end

function path=discover_existing_child(stageRoot)
    path=''; parent=fullfile(stageRoot,'revised_ami_vs_snr_checkpointed_flat');
    if exist(parent,'dir')~=7,return;end
    D=dir(parent); D=D([D.isdir]);
    names={D.name}; keep=~ismember(names,{'.','..'}); D=D(keep);
    if isempty(D),return;end
    if numel(D)>1
        error('sim_i11p1_production_freeze:AmbiguousStageResume', ...
            'More than one child run exists under %s.',parent);
    end
    path=fullfile(parent,D(1).name);
end

function save_manifest_atomic(path,manifest)
    folder=fileparts(path); tmp=[tempname(folder) '.mat'];
    cleanup=onCleanup(@() cleanup_tmp(tmp)); %#ok<NASGU>
    save(tmp,'manifest','-v7');
    [ok,msg]=movefile(tmp,path,'f'); if ~ok,error('sim_i11p1_production_freeze:ManifestSave','%s',msg);end
end
function cleanup_tmp(path),if exist(path,'file')==2,try,delete(path);catch,end,end,end
function token=sanitize_label(v),token=strtrim(char(string(v)));if isempty(token),token='run';end;token=regexprep(token,'[^A-Za-z0-9._-]+','_');end
function tf=text_scalar(v),tf=ischar(v)||(isstring(v)&&isscalar(v));end
function tf=logical_scalar(v),tf=islogical(v)&&isscalar(v);end
function tf=positive_integer(v),tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&isreal(v)&&v>0&&mod(v,1)==0;end
