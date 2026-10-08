function combined = merge_k4_min_gap_results(originalInput,extensionInput,varargin)
%MERGE_K4_MIN_GAP_RESULTS  Merge validated Commit K + K4 spacing data.
%
% The 10 original K physical points and 6 new K4 points form 16 points.
% Inputs can be results structs or paths to their MAT files.
% Existing files are NEVER modified. Derived gain-retention values for
% the six new points are computed using the original d_min=0 baseline.
%
%   C4=merge_k4_min_gap_results(originalMat,R4);
%
% By default output is under <K4 run>/combined and includes one combined
% figure plus three result tables and MAT results.
%
% No SA/AMI evaluation or global optimum claim is made.

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'OutputDirectory','',@text_scalar);
    addParameter(p,'MakeFigure',true,@logical_scalar);
    addParameter(p,'WriteCSV',true,@logical_scalar);
    parse(p,varargin{:}); o=p.Results;

    old=load_result(originalInput);
    ext=load_result(extensionInput);
    verify_compatibility(old,ext);

    if strlength(string(o.OutputDirectory))>0
        outDir=char(string(o.OutputDirectory));
    elseif isfield(ext,'outputDirectory') && ~isempty(ext.outputDirectory)
        outDir=fullfile(ext.outputDirectory,'combined');
    else
        error('merge_k4_min_gap_results:OutputDirectory', ...
            'Provide OutputDirectory when extension results do not include an output directory.');
    end
    ensure_revision2_path(outDir);
    if exist(outDir,'dir')~=7
        [ok,msg]=mkdir(outDir);
        if ~ok,error('merge_k4_min_gap_results:mkdir','%s',msg);end
    end

    restartSummary=sortrows([old.restartSummary;ext.restartSummary], ...
        {'M','MinGap','Replicate','Restart'});
    replicateSummary=sortrows([old.replicateSummary;ext.replicateSummary], ...
        {'M','MinGap','Replicate'});
    physicalSummary=sortrows([old.physicalSummary;ext.physicalSummary], ...
        {'M','MinGap'});

    % The extension has no d_min=0 rows and cannot independently compute
    % retention. Use the original d_min=0 gain separately for each M.
    for m=old.config.MVec
        idx0=old.physicalSummary.M==m & ...
            abs(old.physicalSummary.MinGap)<1e-12;
        if nnz(idx0)~=1
            error('merge_k4_min_gap_results:Reference', ...
                'Expected exactly one original d_min=0 point for M=%d.',m);
        end
        g0=old.physicalSummary.MeanGainOverPAM(idx0);
        idxNew=physicalSummary.M==m & ...
            ~ismembertol(physicalSummary.MinGap,old.config.MinGapVec,1e-10);
        physicalSummary.GainRetentionVsGap0(idxNew)= ...
            physicalSummary.MeanGainOverPAM(idxNew)/g0;
    end

    validate_merged_tables(old,ext,restartSummary,replicateSummary,physicalSummary);

    C=old.config;
    C.version='K4-merged-v1';
    C.description='Combined original Commit K and six feasible Commit K4 extension points';
    C.MinGapVec=unique([old.config.MinGapVec ext.config.MinGapVec]);
    C.MinGapByM=cell(size(C.MVec));
    for i=1:numel(C.MVec)
        m=C.MVec(i);
        C.MinGapByM{i}=physicalSummary.MinGap(physicalSummary.M==m).';
    end
    C.nPhysicalSensitivityPoints=height(physicalSummary);
    C.nReplicateCases=height(replicateSummary);
    C.nRestartTasks=height(restartSummary);
    C.totalSAIterations=old.config.totalSAIterations+ext.config.totalSAIterations;
    C.nominalFastObjectiveEvaluations= ...
        old.config.nominalFastObjectiveEvaluations+ext.config.nominalFastObjectiveEvaluations;

    combined=struct();
    combined.version='K4-merged-v1';
    combined.config=C;
    combined.baselines=old.baselines;
    combined.restartSummary=restartSummary;
    combined.replicateSummary=replicateSummary;
    combined.physicalSummary=physicalSummary;
    combined.outputDirectory=outDir;
    combined.sources=struct('original',source_label(originalInput), ...
        'extension',source_label(extensionInput));
    combined.metrics=struct( ...
        'nPhysicalSensitivityPoints',height(physicalSummary), ...
        'nReplicateCases',height(replicateSummary), ...
        'nRestartTasks',height(restartSummary), ...
        'newRestartTasks',height(ext.restartSummary), ...
        'maxAbsFastValidatorGap',max(abs(restartSummary.FastValidatorGap)), ...
        'selectionChanges',sum(replicateSummary.SelectionChanged), ...
        'maxReplicateRange',max(physicalSummary.ReplicateRange), ...
        'originalScheduler',old.metrics.scheduler, ...
        'extensionScheduler',ext.metrics.scheduler);
    combined.globalOptimumClaim=false;

    if o.WriteCSV
        writetable(restartSummary,fullfile(outDir,'K4_combined_restart_results.csv'));
        writetable(replicateSummary,fullfile(outDir,'K4_combined_replicate_summary.csv'));
        writetable(physicalSummary,fullfile(outDir,'K4_combined_physical_summary.csv'));
    end

    if o.MakeFigure
        combined.figurePaths=plot_k4_combined_min_gap(combined, ...
            'OutputDirectory',outDir);
    else
        combined.figurePaths=struct();
    end

    result=combined; %#ok<NASGU>
    save(fullfile(outDir,'K4_combined_results.mat'),'result','-v7.3');

    fprintf('\n============================================================\n');
    fprintf('COMMIT K4 - MERGED RESULTS COMPLETE\n');
    fprintf('Original points: %d | new points: %d | combined points: %d\n', ...
        height(old.physicalSummary),height(ext.physicalSummary), ...
        height(physicalSummary));
    fprintf('Replicates: %d | restart tasks: %d | new restarts: %d\n', ...
        height(replicateSummary),height(restartSummary),height(ext.restartSummary));
    fprintf('Max fast-validator gap: %.3e bit/symbol\n', ...
        combined.metrics.maxAbsFastValidatorGap);
    fprintf('Output directory: %s\n',outDir);
    fprintf('M/d_min/validated AMI (two-replicate mean):\n');
    disp(physicalSummary(:,{'M','MinGap','Rep1ValidatedAMI', ...
        'Rep2ValidatedAMI','MeanValidatedAMI','StdAcrossReplicates', ...
        'MeanGainOverPAM','GainRetentionVsGap0'}));
    fprintf('Original Commit-K files are unchanged.\n');
    fprintf('No global-optimum claim is implied.\n');
    fprintf('============================================================\n');
end


function verify_compatibility(old,ext)
    for v={'config','baselines','restartSummary','replicateSummary','physicalSummary','metrics'}
        f=v{1};
        if ~isfield(old,f)||~isfield(ext,f)
            error('merge_k4_min_gap_results:MissingField','Both results require %s.',f);
        end
    end

    if ~strcmp(old.version,'K-v1')||~strcmp(ext.version,'K4-v1')
        error('merge_k4_min_gap_results:Versions', ...
            'Require frozen original K-v1 and extension K4-v1 results.');
    end

    common={'MVec','SNRdB','SigmaR2','Replicates','NStarts','MaxIter', ...
        'P_avg','pinZero','baseSeed','FastFadingMethod','LoghDtPolicy', ...
        'loghDt','loghDtRefined','loghSpanSigma','loghYBlockSize', ...
        'YGridMode','YGridScale','historyEvery','T0','Tf', ...
        'BaseStd0','ItersPerTemp','SelectedDtByM'};
    for k=1:numel(common)
        f=common{k};
        if ~isfield(old.config,f)||~isfield(ext.config,f)|| ...
                ~isequaln(old.config.(f),ext.config.(f))
            error('merge_k4_min_gap_results:ConfigMismatch', ...
                'Scientific configuration mismatch: %s.',f);
        end
    end

    originalExpected=[0 0.005 0.01 0.02 0.05];
    extensionExpected=[0.03 0.04 0.07 0.10];
    if ~isequal(old.config.MinGapVec,originalExpected)|| ...
            ~isequal(ext.config.MinGapVec,extensionExpected)|| ...
            ~isequal(ext.config.MinGapByM, ...
            {[0.03 0.04 0.07 0.10],[0.03 0.04]})
        error('merge_k4_min_gap_results:GapProtocol','Unexpected gap grids.');
    end

    for i=1:numel(old.config.MVec)
        a=old.baselines(i); b=ext.baselines(i);
        if a.M~=b.M || abs(a.amiValidated-b.amiValidated)>1e-9 || ...
                abs(a.amiFast-b.amiFast)>1e-9 || ...
                max(abs(a.x(:)-b.x(:)))>1e-12
            error('merge_k4_min_gap_results:PAMBaseline', ...
                'Uniform-PAM baselines disagree for M=%d.',a.M);
        end
    end

    for m=old.config.MVec
        for rep=old.config.Replicates
            x=unique(old.casePlan.Seed(old.casePlan.M==m & ...
                old.casePlan.Replicate==rep));
            y=unique(ext.casePlan.Seed(ext.casePlan.M==m & ...
                ext.casePlan.Replicate==rep));
            if numel(x)~=1||numel(y)~=1||x~=y
                error('merge_k4_min_gap_results:SeedMismatch', ...
                    'Paired seeds disagree for M=%d replicate=%d.',m,rep);
            end
        end
    end
end


function validate_merged_tables(old,ext,T,R,P)
    if height(old.restartSummary)~=160||height(ext.restartSummary)~=96|| ...
            height(T)~=256||height(R)~=32||height(P)~=16
        error('merge_k4_min_gap_results:Counts', ...
            'Expected 160+96=256 restarts, 32 replicates, 16 physical points.');
    end

    if ~isequaln(sortrows(T(ismembertol(T.MinGap,old.config.MinGapVec,1e-10),:), ...
            {'M','MinGap','Replicate','Restart'}), ...
            sortrows(old.restartSummary,{'M','MinGap','Replicate','Restart'}))|| ...
       ~isequaln(sortrows(R(ismembertol(R.MinGap,old.config.MinGapVec,1e-10),:), ...
            {'M','MinGap','Replicate'}), ...
            sortrows(old.replicateSummary,{'M','MinGap','Replicate'}))|| ...
       ~isequaln(sortrows(P(ismembertol(P.MinGap,old.config.MinGapVec,1e-10),:), ...
            {'M','MinGap'}),sortrows(old.physicalSummary,{'M','MinGap'}))
        error('merge_k4_min_gap_results:OriginalChanged', ...
            'Original K table rows changed during merge.');
    end

    keys=unique(T(:,{'M','MinGap','Replicate','Restart'}),'rows','sorted');
    if height(keys)~=256
        error('merge_k4_min_gap_results:DuplicateRestart','Duplicate restart key.');
    end
    keys=unique(R(:,{'M','MinGap','Replicate'}),'rows','sorted');
    if height(keys)~=32
        error('merge_k4_min_gap_results:DuplicateReplicate','Duplicate replicate key.');
    end
    keys=unique(P(:,{'M','MinGap'}),'rows','sorted');
    if height(keys)~=16
        error('merge_k4_min_gap_results:DuplicatePhysical','Duplicate physical point.');
    end

    for i=1:height(P)
        m=P.M(i);
        d=P.MinGap(i);
        if d>2*old.config.P_avg/(m-1)+1e-12
            error('merge_k4_min_gap_results:Infeasible', ...
                'M=%d, d_min=%g violates fixed mean optical power.',m,d);
        end
        row=R(R.M==m & abs(R.MinGap-d)<1e-12,:);
        if height(row)~=2 || ~isequal(sort(row.Replicate).',[1 2])
            error('merge_k4_min_gap_results:ReplicateCount', ...
                'M=%d, d_min=%g missing independent replicate.',m,d);
        end
        if abs(mean(row.WinnerValidatedAMI)-P.MeanValidatedAMI(i))>1e-9 || ...
                abs(std(row.WinnerValidatedAMI,0)-P.StdAcrossReplicates(i))>1e-9
            error('merge_k4_min_gap_results:SummaryMismatch', ...
                'AMI/statistics mismatch for M=%d, d_min=%g.',m,d);
        end
        for j=1:height(row)
            rr=T(T.M==m & abs(T.MinGap-d)<1e-12 & ...
                T.Replicate==row.Replicate(j),:);
            if height(rr)~=8 || ~any(rr.Restart==row.WinnerRestart(j)) || ...
                    abs(max(rr.ValidatedAMI)-row.WinnerValidatedAMI(j))>1e-9
                error('merge_k4_min_gap_results:Restarts', ...
                    'Restart/winner inconsistency at M=%d d_min=%g.',m,d);
            end
            x=double(row.WinnerConstellation{j}(:));
            if numel(x)~=m || any(diff(x)<d-1e-8) || ...
                    abs(mean(x)-old.config.P_avg)>1e-8 || abs(x(1))>1e-8
                error('merge_k4_min_gap_results:Geometry', ...
                    'Selected constellation violates constraints at M=%d d_min=%g.',m,d);
            end
        end
    end
end


function result=load_result(v)
    if isstruct(v)
        result=v;return;
    end
    if ~(ischar(v)||(isstring(v)&&isscalar(v)))
        error('merge_k4_min_gap_results:Input', ...
            'Input must be a result struct or a MAT-file path.');
    end
    path=char(string(v));
    if exist(path,'file')~=2
        error('merge_k4_min_gap_results:MissingFile','Missing MAT file: %s',path);
    end
    S=load(path,'result');
    if ~isfield(S,'result')
        error('merge_k4_min_gap_results:Variable', ...
            'MAT-file must contain result struct.');
    end
    result=S.result;
end


function str=source_label(v)
    if isstruct(v)
        if isfield(v,'outputDirectory')
            str=char(string(v.outputDirectory));
        else
            str='in-memory result';
        end
    else
        str=char(string(v));
    end
end


function ensure_revision2_path(path)
    root=char(java.io.File(fso_result_utils.revision2_results_root()).getCanonicalPath());
    target=char(java.io.File(path).getCanonicalPath());
    if ~(strcmp(target,root)||startsWith(target,[root filesep]))
        error('merge_k4_min_gap_results:ResultsRoot', ...
            'Combined outputs must be under results/ieee_revision2.');
    end
end


function tf=text_scalar(v),tf=ischar(v)||(isstring(v)&&isscalar(v));end
function tf=logical_scalar(v),tf=islogical(v)&&isscalar(v);end
