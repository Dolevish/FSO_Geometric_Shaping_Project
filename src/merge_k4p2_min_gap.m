function C42 = merge_k4p2_min_gap(k4CombinedInput,k4p2Input,varargin)
%MERGE_K4P2_MIN_GAP  Replace d_min=0 rows with new d_min=0.001 rows.
%
% Immutable source K4 combined results contain 256 restarts/16 points:
% original 0/.005/.01/.02/.03/.04/.05, plus .07/.10 for M16.
% K4.2 adds exactly 32 new restarts for M16/M32 at .001 (2 replicates,
% 8 starts each). Original .005 is deliberately REUSED, never recomputed.
% Final derived dataset retains 256 restarts, 32 reps, 16 physical points.
%
%   C42=merge_k4p2_min_gap('/path/K4_combined_results.mat',R42);
%
% Existing source MAT files are NEVER modified. Output tables explicitly
% use GainRetentionVsGap0001 (reference=.001), NOT stale d_min=0 gain.
%
% To apply the exact user FIG styling, pass 'TemplateFigPath','/path/*.fig'.

    p=inputParser;p.FunctionName=mfilename;
    addParameter(p,'OutputDirectory','',@is_text);
    addParameter(p,'TemplateFigPath','',@is_text);
    addParameter(p,'MakeFigure',true,@is_bool);
    addParameter(p,'WriteCSV',true,@is_bool);
    parse(p,varargin{:}); o=p.Results;

    old=read_result(k4CombinedInput);
    fresh=read_result(k4p2Input);
    verify_scientific_sources(old,fresh);

    if strlength(string(o.OutputDirectory))>0
        outDir=char(string(o.OutputDirectory));
    elseif isfield(fresh,'outputDirectory') && ...
            ~isempty(fresh.outputDirectory)
        outDir=fullfile(char(fresh.outputDirectory),'combined');
    else
        error('K4p2:MissingOutput','Specify OutputDirectory.');
    end
    root=char(java.io.File(fso_result_utils.revision2_results_root()).getCanonicalPath());
    target=char(java.io.File(outDir).getCanonicalPath());
    if ~(strcmp(root,target)||startsWith(target,[root filesep]))
        error('K4p2:OutputRoot','Outputs must be under results/ieee_revision2.');
    end

    % Existing results with zero removed. The ORIGINAL source remains intact.
    nRestart=old.restartSummary(abs(old.restartSummary.MinGap)>1e-12,:);
    nRep=old.replicateSummary(abs(old.replicateSummary.MinGap)>1e-12,:);
    nPhys=old.physicalSummary(abs(old.physicalSummary.MinGap)>1e-12,:);
    if height(nRestart)~=224 || height(nRep)~=28 || height(nPhys)~=14
        error('K4p2:HistoricalCounts','Expected to remove exactly 32/4/2 zero-gap rows.');
    end
    ensure_no_overlap(nRestart,fresh.restartSummary, ...
        {'M','MinGap','Replicate','Restart'});
    ensure_no_overlap(nRep,fresh.replicateSummary, ...
        {'M','MinGap','Replicate'});
    ensure_no_overlap(nPhys,fresh.physicalSummary,{'M','MinGap'});
    restartSummary=sortrows([nRestart;fresh.restartSummary], ...
        {'M','MinGap','Replicate','Restart'});
    replicateSummary=sortrows([nRep;fresh.replicateSummary], ...
        {'M','MinGap','Replicate'});
    physicalSummary=sortrows([nPhys;fresh.physicalSummary],{'M','MinGap'});

    % Recompute gain retention using the smallest INCLUDED gap .001.
    % The column GainRetentionVsGap0 becomes scientifically invalid once
    % zero-gap rows are omitted; delete it and replace with honest metadata.
    if ~ismember('GainRetentionVsGap0',physicalSummary.Properties.VariableNames)
        error('K4p2:Schema','Original gain-retention column missing.');
    end
    physicalSummary.GainRetentionVsGap0=[];
    physicalSummary.GainRetentionVsGap0001=zeros(height(physicalSummary),1);
    for m=[16 32]
        ix=physicalSummary.M==m;
        iRef=ix & abs(physicalSummary.MinGap-.001)<1e-12;
        if nnz(iRef)~=1,error('K4p2:Reference','Missing .001 M=%d.',m);end
        ref=physicalSummary.MeanGainOverPAM(iRef);
        if ~isfinite(ref)||ref<=0
            error('K4p2:ReferenceGain','Non-positive .001 reference gain.');
        end
        physicalSummary.GainRetentionVsGap0001(ix)= ...
            physicalSummary.MeanGainOverPAM(ix)/ref;
    end

    validate_final(old,fresh,nRestart,nRep,nPhys, ...
        restartSummary,replicateSummary,physicalSummary);

    cfg=old.config;
    cfg.version='K4p2-merged-v1';
    cfg.description='K4+K4.2: original d_min=0 excluded and 0.001 added';
    cfg.MinGapVec=unique(physicalSummary.MinGap).';
    cfg.MinGapByM=cell(size(cfg.MVec));
    for i=1:numel(cfg.MVec)
        cfg.MinGapByM{i}=physicalSummary.MinGap( ...
            physicalSummary.M==cfg.MVec(i)).';
    end
    cfg.nPhysicalSensitivityPoints=height(physicalSummary);
    cfg.nReplicateCases=height(replicateSummary);
    cfg.nRestartTasks=height(restartSummary);
    cfg.totalSAIterations=height(restartSummary)*cfg.MaxIter;
    cfg.nominalFastObjectiveEvaluations= ...
        height(restartSummary)*(cfg.MaxIter+1);

    C42=struct();
    C42.version='K4p2-merged-v1';
    C42.config=cfg;
    C42.baselines=old.baselines;
    C42.restartSummary=restartSummary;
    C42.replicateSummary=replicateSummary;
    C42.physicalSummary=physicalSummary;
    C42.sources=struct('K4Combined',source_label(k4CombinedInput), ...
        'K4p2New',source_label(k4p2Input));
    C42.metrics=struct('originalNonzeroRestarts',224, ...
        'originalZeroRemoved',32,'newRestarts',32, ...
        'totalRestartTasks',256,'replicateCases',32,'physicalPoints',16, ...
        'maximumFastValidatorGap',max(abs(restartSummary.FastValidatorGap)), ...
        'winnerSelectionChanges',sum(replicateSummary.SelectionChanged), ...
        'referenceMinGap',.001,'globalOptimumClaim',false);
    C42.outputDirectory=outDir;
    C42.globalOptimumClaim=false;

    if exist(outDir,'dir')~=7
        [ok,msg]=mkdir(outDir);
        if ~ok,error('K4p2:mkdir','%s',msg);end
    end
    if o.WriteCSV
        writetable(restartSummary,fullfile(outDir, ...
            'K4p2_combined_restart_results.csv'));
        writetable(replicateSummary,fullfile(outDir, ...
            'K4p2_combined_replicate_summary.csv'));
        writetable(physicalSummary,fullfile(outDir, ...
            'K4p2_combined_physical_summary.csv'));
    end
    if o.MakeFigure
        C42.figurePaths=plot_k4p2_min_gap(C42, ...
            'OutputDirectory',outDir, ...
            'TemplateFigPath',o.TemplateFigPath);
    else
        C42.figurePaths=struct();
    end
    result=C42; %#ok<NASGU>
    save(fullfile(outDir,'K4p2_combined_results.mat'),'result','-v7.3');

    fprintf('\n============================================================\n');
    fprintf('COMMIT K4.2 - MERGED RESULTS COMPLETE\n');
    fprintf('Source K4: 256 restarts; removed 32 gap=0; added 32 gap=.001\n');
    fprintf('Final: 16 physical points | 32 replicates | 256 restarts\n');
    fprintf('d_min=0.005 retained from original production; NO rerun.\n');
    fprintf('Reference for gain retention is now d_min=0.001.\n');
    fprintf('Fast-validator max abs gap: %.3e bits/symbol\n', ...
        C42.metrics.maximumFastValidatorGap);
    disp(physicalSummary(:,{'M','MinGap','Rep1ValidatedAMI', ...
        'Rep2ValidatedAMI','MeanValidatedAMI','StdAcrossReplicates', ...
        'MeanGainOverPAM','GainRetentionVsGap0001'}));
    fprintf('Output: %s\n',outDir);
    fprintf('Original K and K4 files are unchanged; no global optimum claim.\n');
    fprintf('============================================================\n');
end

function verify_scientific_sources(old,fresh)
    if ~strcmp(old.version,'K4-merged-v1')||~strcmp(fresh.version,'K4p2-v1')
        error('K4p2:Versions','Require K4 merged results and K4p2 extension.');
    end
    fields={'MVec','SNRdB','SigmaR2','Replicates','NStarts','MaxIter', ...
        'P_avg','pinZero','baseSeed','FastFadingMethod','LoghDtPolicy', ...
        'loghDt','loghDtRefined','loghSpanSigma','loghYBlockSize', ...
        'YGridMode','YGridScale','historyEvery','T0','Tf','BaseStd0', ...
        'ItersPerTemp','SelectedDtByM'};
    for k=1:numel(fields)
        f=fields{k};
        if ~isfield(old.config,f)||~isfield(fresh.config,f)|| ...
                ~isequaln(old.config.(f),fresh.config.(f))
            error('K4p2:FrozenMismatch','Mismatch in scientific field %s.',f);
        end
    end
    if ~isequal(old.config.MVec,[16 32])||old.config.SNRdB~=20|| ...
            old.config.SigmaR2~=.2||old.config.P_avg~=1|| ...
            ~old.config.pinZero||~isequal(old.config.Replicates,[1 2])|| ...
            old.config.NStarts~=8||old.config.MaxIter~=8000|| ...
            ~isequal(fresh.config.MinGapVec,.001)|| ...
            ~isequal(fresh.config.MinGapByM,{.001,.001}) || ...
            any(abs(old.config.SelectedDtByM-.01)>1e-12)
        error('K4p2:FrozenProtocol','K4.2 must use the exact frozen K policy.');
    end
    if height(old.restartSummary)~=256|| ...
            height(old.replicateSummary)~=32|| ...
            height(old.physicalSummary)~=16|| ...
            height(fresh.restartSummary)~=32|| ...
            height(fresh.replicateSummary)~=4|| ...
            height(fresh.physicalSummary)~=2
        error('K4p2:SourcesIncomplete','Expected complete K4 and K4.2 sources.');
    end
    for m=[16 32]
        a=old.baselines([old.baselines.M]==m);
        b=fresh.baselines([fresh.baselines.M]==m);
        if numel(a)~=1||numel(b)~=1 || ...
                abs(a.amiValidated-b.amiValidated)>1e-9 || ...
                abs(a.amiFast-b.amiFast)>1e-9 || ...
                max(abs(a.x(:)-b.x(:)))>1e-12
            error('K4p2:PAMMismatch','PAM baseline mismatch at M=%d.',m);
        end
        for rep=[1 2]
            oldSeed=unique(old.replicateSummary.Seed( ...
                old.replicateSummary.M==m & ...
                old.replicateSummary.Replicate==rep));
            newSeed=unique(fresh.replicateSummary.Seed( ...
                fresh.replicateSummary.M==m & ...
                fresh.replicateSummary.Replicate==rep));
            if numel(oldSeed)~=1 || numel(newSeed)~=1 || oldSeed~=newSeed
                error('K4p2:Seed','Paired seed mismatch M=%d rep=%d.',m,rep);
            end
        end
    end
end

function ensure_no_overlap(A,B,fields)
    a=A(:,fields);b=B(:,fields);
    both=[a;b];
    if height(unique(both,'rows'))~=height(both)
        error('K4p2:DuplicateCase','Duplicate task/case in merged files.');
    end
end

function validate_final(old,fresh,unchangedT,unchangedR,unchangedP,T,R,P)
    if height(T)~=256||height(R)~=32||height(P)~=16 || ...
            height(unique(T(:,{'M','MinGap','Replicate','Restart'}),'rows'))~=256 || ...
            height(unique(R(:,{'M','MinGap','Replicate'}),'rows'))~=32 || ...
            height(unique(P(:,{'M','MinGap'}),'rows'))~=16 || ...
            any(abs(T.MinGap)<1e-12)||any(abs(R.MinGap)<1e-12)|| ...
            any(abs(P.MinGap)<1e-12)
        error('K4p2:UniqueCounts','Unexpected final counts/duplicates/zero gaps.');
    end
    for m=[16 32]
        x=P.MinGap(P.M==m);
        expectation=[.001 .005 .01 .02 .03 .04 .05];
        if m==16,expectation=[expectation .07 .10];end
        if numel(x)~=numel(expectation) || ...
                any(abs(x(:)-expectation(:))>1e-12)
            error('K4p2:MissingPoint','Missing/extra distance point M=%d.',m);
        end
    end
    if ~isequaln(sortrows(T(abs(T.MinGap-.001)>1e-12,:), ...
            {'M','MinGap','Replicate','Restart'}),unchangedT)|| ...
       ~isequaln(sortrows(R(abs(R.MinGap-.001)>1e-12,:), ...
            {'M','MinGap','Replicate'}),unchangedR)
        error('K4p2:HistoricalChange','Historical restart/replicate rows changed.');
    end
    historicalP=P(abs(P.MinGap-.001)>1e-12,:);
    historicalP.GainRetentionVsGap0001=[];
    unchangedP.GainRetentionVsGap0=[];
    if ~isequaln(historicalP,unchangedP)
        error('K4p2:HistoricalPhysical','Historical AMI physical values changed.');
    end
    for k=1:height(P)
        m=P.M(k);d=P.MinGap(k);
        if d>2*old.config.P_avg/(m-1)+1e-12
            error('K4p2:Feasibility','Infeasible d_min M=%d.',m);
        end
        r=R(R.M==m & abs(R.MinGap-d)<1e-12,:);
        if height(r)~=2||~isequal(sort(r.Replicate).',[1 2])
            error('K4p2:ReplicateCount','Expected exactly 2 replicates.');
        end
        if abs(mean(r.WinnerValidatedAMI)-P.MeanValidatedAMI(k))>1e-9 || ...
                abs(std(r.WinnerValidatedAMI,0)-P.StdAcrossReplicates(k))>1e-9
            error('K4p2:AMIStats','Replicate AMI aggregation mismatch.');
        end
        for j=1:2
            rr=T(T.M==m & abs(T.MinGap-d)<1e-12 & ...
                T.Replicate==r.Replicate(j),:);
            if height(rr)~=8 || ...
                    abs(max(rr.ValidatedAMI)-r.WinnerValidatedAMI(j))>1e-9 || ...
                    ~any(rr.Restart==r.WinnerRestart(j))
                error('K4p2:RestartWinner','Missing validated winning restart.');
            end
            x=r.WinnerConstellation{j}(:);
            if numel(x)~=m || any(diff(x)<d-1e-7) || ...
                    abs(mean(x)-old.config.P_avg)>1e-8 || abs(x(1))>1e-8
                error('K4p2:Constellation','Winner violates constellation constraints.');
            end
        end
    end
    if ~isequaln(P.GainRetentionVsGap0001(P.MinGap==.001),ones(2,1))
        error('K4p2:GainReference','Reference .001 gain must be 100%%.');
    end
end

function value=read_result(v)
    if isstruct(v),value=v;return;end
    if ~(ischar(v)||isstring(v))
        error('K4p2:Input','Result must be a struct or a MAT filename.');
    end
    S=load(char(string(v)),'result');
    if ~isfield(S,'result')
        error('K4p2:MATField','Source MAT has no result variable.');
    end
    value=S.result;
end
function s=source_label(v)
    if isstruct(v)
        if isfield(v,'outputDirectory'),s=char(string(v.outputDirectory));
        else,s='in-memory input';end
    else
        s=char(string(v));
    end
end
function t=is_text(v),t=ischar(v)||(isstring(v)&&isscalar(v));end
function t=is_bool(v),t=islogical(v)&&isscalar(v);end
