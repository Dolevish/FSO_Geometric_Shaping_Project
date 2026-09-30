function d = diagnose_i10p2p1_runtime_tuning_sweep(varargin)
%DIAGNOSE_I10P2P1_RUNTIME_TUNING_SWEEP  Production-grid runtime tuning sweep.
%
% Commit I.10.2.1 broadens the I.10.2 kernel benchmark across representative
% production conditions before any runtime-tuning policy is frozen.
%
% The default case matrix is [M SNR_dB sigma_X_sq dt].  dt is explicit so
% the sweep can include the anticipated refined boundary case
% (M=16,SNR=30,sigma_X_sq=0.2,dt=0.005) without changing select_logh_dt().
% A 3-column case matrix [M SNR sigma] is also accepted; then dt is selected
% by the current select_logh_dt() policy.
%
% Two geometries are tested by default:
%   uniform - deterministic uniform positive M-PAM baseline
%   wide    - upper-eighth-style stress geometry with the same mean power
%
% Every scenario compares:
%   reference: exact I.10.1 cached kernel, fixed block=512, batch=1
%   auto     : current I.10.2 runtime policy
%   tuned    : explicit block-size x symbol-batch combinations using I.10.2
%
% No SA is run. No scientific setting or production selector is changed.

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'Cases',default_cases(),@cases_valid);
    addParameter(p,'GeometryModes',["uniform","wide"],@geometry_modes_valid);
    addParameter(p,'BlockSizes',[256 512 1024 2048],@positive_integer_vector);
    addParameter(p,'BatchSizes',[1 2 4],@positive_integer_vector);
    addParameter(p,'nWarmup',1,@nonnegative_integer);
    addParameter(p,'nRepeats',3,@positive_integer);
    addParameter(p,'AdaptiveScale',1.1,@(v) isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>=1);
    addParameter(p,'minGap',0.01,@nonnegative_scalar);
    addParameter(p,'AbsAgreementTolerance',5e-12,@positive_scalar);
    addParameter(p,'NoRegressionThreshold',0.98,@positive_scalar);
    addParameter(p,'WriteCSV',true,@(v) islogical(v)&&isscalar(v));
    addParameter(p,'RunLabel','',@(v) ischar(v)||(isstring(v)&&isscalar(v)));
    addParameter(p,'ResultsRoot',fso_result_utils.default_results_root(),@(v) ischar(v)||isstring(v));
    parse(p,varargin{:}); o=p.Results;

    C=double(o.Cases);
    geometries=unique(lower(string(o.GeometryModes(:).')),'stable');
    blocks=unique(double(o.BlockSizes(:).'),'stable');
    batches=unique(double(o.BatchSizes(:).'),'stable');

    parent=fullfile(char(o.ResultsRoot),'i10p2p1_runtime_tuning');
    label=sanitize_label(o.RunLabel);
    utc=char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmssSSS''Z'''));
    [outDir,runId]=make_unique_dir(parent,sprintf('%s_%s',label,utc));

    fprintf('\n============================================================\n');
    fprintf('Commit I.10.2.1 - production-grid runtime tuning sweep\n');
    fprintf('Cases: %d | geometries=%s\n',size(C,1),strjoin(cellstr(geometries),','));
    fprintf('Blocks=%s | batches=%s | repeats=%d | warmup=%d\n', ...
        mat2str(blocks),mat2str(batches),o.nRepeats,o.nWarmup);
    fprintf('No-regression threshold: %.3fx\n',o.NoRegressionThreshold);
    fprintf('Results: %s\n',outDir);
    fprintf('============================================================\n\n');

    rows={}; q=0; scenario=0;
    for k=1:size(C,1)
        M=C(k,1); snr=C(k,2); sig=C(k,3);
        if size(C,2)==4
            dt=C(k,4); dtSource="explicit";
        else
            [dt,~]=select_logh_dt(M,snr,sig); dtSource="selector";
        end

        cfg=build_fso_config(M,1,snr,sig,'minGap',o.minGap,'pinZero',true, ...
            'saMaxIter',1,'saNStarts',1,'YGridMode','candidate-adaptive', ...
            'YGridScale',o.AdaptiveScale,'loghDt',dt,'loghSpanSigma',8,'loghYBlockSize',512);

        for g=1:numel(geometries)
            geometry=geometries(g); scenario=scenario+1;
            x=make_geometry(geometry,M,1,o.minGap);
            x=x(:); AMI_functions.assert_constellation_feasible(x,cfg,1e-10);

            % Exact I.10.1 cached implementation reference.
            ref=cfg; ref.loghRuntimeKernel='i10p1'; ref.loghRuntimeBlockPolicy='fixed';
            ref.loghSymbolBatchSize=1; ref.loghZeroFastPath=false;
            warm_eval(x,ref,dt,o.AdaptiveScale,512,o.nWarmup);
            [refAMI,refMed,refMean,refMin]=timed_eval(x,ref,dt,o.AdaptiveScale,512,o.nRepeats);

            % Current automatic I.10.2 policy.
            auto=cfg;
            [~,dAuto]=AMI_noCSI_fast_logh_adaptive_grid(x,auto.px,auto,dt,8,o.AdaptiveScale,512);
            warm_eval(x,auto,dt,o.AdaptiveScale,512,o.nWarmup);
            [v,tmed,tmean,tmin]=timed_eval(x,auto,dt,o.AdaptiveScale,512,o.nRepeats);
            assert_agreement(v,refAMI,o.AbsAgreementTolerance,M,snr,sig,geometry,'auto');
            q=q+1; rows{q,1}=make_row("auto",scenario,k,geometry,M,sig,snr,dt,dtSource,max(x), ...
                512,dAuto.yBlockSizeEffective,dAuto.symbolBatchSize,refAMI,refMed,refMean,refMin, ...
                v,tmed,tmean,tmin,abs(v-refAMI)); %#ok<AGROW>

            % Explicit implementation-only sweep.
            for b=blocks
                for sb=batches
                    if sb>M, continue; end
                    tuned=cfg; tuned.loghRuntimeKernel='i10p2'; tuned.loghRuntimeBlockPolicy='fixed';
                    tuned.loghSymbolBatchSize=sb; tuned.loghZeroFastPath=true;
                    [~,dTune]=AMI_noCSI_fast_logh_adaptive_grid(x,tuned.px,tuned,dt,8,o.AdaptiveScale,b);
                    warm_eval(x,tuned,dt,o.AdaptiveScale,b,o.nWarmup);
                    [v,tmed,tmean,tmin]=timed_eval(x,tuned,dt,o.AdaptiveScale,b,o.nRepeats);
                    labelTune=sprintf('b%d-s%d',b,sb);
                    assert_agreement(v,refAMI,o.AbsAgreementTolerance,M,snr,sig,geometry,labelTune);
                    q=q+1; rows{q,1}=make_row("tuned",scenario,k,geometry,M,sig,snr,dt,dtSource,max(x), ...
                        b,dTune.yBlockSizeEffective,dTune.symbolBatchSize,refAMI,refMed,refMean,refMin, ...
                        v,tmed,tmean,tmin,abs(v-refAMI)); %#ok<AGROW>
                end
            end
        end
        fprintf('  completed %d/%d physical cases\n',k,size(C,1));
    end

    T=struct2table(vertcat(rows{:}));
    tuned=T(T.Mode=="tuned",:);
    auto=T(T.Mode=="auto",:);

    aggregate=aggregate_tuned_by_M(tuned);
    autoAgg=aggregate_auto_by_M(auto);
    [recommended,bestMedian]=recommend_by_M(aggregate,o.NoRegressionThreshold);

    summary=struct();
    summary.nPhysicalCases=size(C,1);
    summary.nScenarios=scenario;
    summary.nRows=height(T);
    summary.maxAbsDifference=max(T.AbsDifference,[],'omitnan');
    summary.noRegressionThreshold=o.NoRegressionThreshold;
    summary.autoAggregateByM=autoAgg;
    summary.tuningAggregate=aggregate;
    summary.recommendedByM=recommended;
    summary.bestMedianByM=bestMedian;

    d=struct();
    d.run=struct('id',runId,'outputDirectory',outDir);
    d.settings=struct('Cases',C,'GeometryModes',geometries,'BlockSizes',blocks, ...
        'BatchSizes',batches,'nWarmup',o.nWarmup,'nRepeats',o.nRepeats, ...
        'AdaptiveScale',o.AdaptiveScale,'minGap',o.minGap, ...
        'AbsAgreementTolerance',o.AbsAgreementTolerance, ...
        'NoRegressionThreshold',o.NoRegressionThreshold);
    d.results=T; d.summary=summary;

    save(fullfile(outDir,'i10p2p1RuntimeTuning.mat'),'d','-v7.3');
    if o.WriteCSV
        writetable(T,fullfile(outDir,'i10p2p1RuntimeTuning_all.csv'));
        writetable(aggregate,fullfile(outDir,'i10p2p1RuntimeTuning_aggregate.csv'));
        writetable(recommended,fullfile(outDir,'i10p2p1RuntimeTuning_recommended.csv'));
    end

    fprintf('\nI.10.2.1 RUNTIME TUNING SUMMARY\n');
    fprintf('physical cases / scenarios / rows : %d / %d / %d\n', ...
        summary.nPhysicalCases,summary.nScenarios,summary.nRows);
    fprintf('max numerical |candidate-I10.1|   : %.3e\n',summary.maxAbsDifference);
    fprintf('\nCurrent auto policy aggregate by M:\n');
    disp(autoAgg);
    fprintf('\nRobust recommendation by M (measured only; selector NOT changed):\n');
    disp(recommended);
    fprintf('\nHighest-median combination by M:\n');
    disp(bestMedian);
end


function C=default_cases()
% Representative low/medium/high production points. The boundary M16 case
% is explicitly evaluated at dt=0.005 because I.10 identified it for likely
% refinement, while leaving the production selector untouched in this commit.
    C=[ ...
         4 10 0.1 0.010; ...
         4 20 0.2 0.010; ...
         4 30 0.3 0.010; ...
         8 10 0.1 0.010; ...
         8 20 0.2 0.010; ...
         8 30 0.3 0.010; ...
        16 10 0.1 0.010; ...
        16 20 0.2 0.010; ...
        16 30 0.1 0.010; ...
        16 30 0.2 0.005; ...
        16 30 0.3 0.005; ...
        32 10 0.1 0.010; ...
        32 20 0.2 0.010; ...
        32 30 0.1 0.005; ...
        32 30 0.2 0.005; ...
        32 30 0.3 0.005];
end


function r=make_row(mode,scenario,caseIndex,geometry,M,sig,snr,dt,dtSource,xMax, ...
    blockReq,blockEff,batch,refAMI,refMed,refMean,refMin,v,tmed,tmean,tmin,err)
    r=struct();
    r.Mode=string(mode); r.Scenario=scenario; r.CaseIndex=caseIndex;
    r.Geometry=string(geometry); r.M=M; r.SigmaX2=sig; r.SNRdB=snr; r.Dt=dt;
    r.DtSource=string(dtSource); r.XMax=xMax;
    r.BlockRequested=blockReq; r.BlockEffective=blockEff; r.BatchSize=batch;
    r.ReferenceAMI=refAMI; r.ReferenceMedianSeconds=refMed;
    r.ReferenceMeanSeconds=refMean; r.ReferenceMinSeconds=refMin;
    r.OptimizedAMI=v; r.MedianOptimizedSeconds=tmed;
    r.MeanOptimizedSeconds=tmean; r.MinOptimizedSeconds=tmin;
    r.ReferenceOverOptimizedSpeedup=refMed/tmed; r.AbsDifference=err;
end


function A=aggregate_tuned_by_M(T)
    Ms=unique(T.M(:).'); rows={}; q=0;
    for M=Ms
        TM=T(T.M==M,:);
        blocks=unique(TM.BlockEffective(:).'); batches=unique(TM.BatchSize(:).');
        for b=blocks
            for sb=batches
                R=TM(TM.BlockEffective==b & TM.BatchSize==sb,:);
                if isempty(R), continue; end
                s=R.ReferenceOverOptimizedSpeedup;
                q=q+1; rows{q,1}=struct( ...
                    'M',M,'BlockSize',b,'BatchSize',sb,'NScenarios',height(R), ...
                    'MedianSpeedup',median(s,'omitnan'),'MeanSpeedup',mean(s,'omitnan'), ...
                    'MinSpeedup',min(s,[],'omitnan'),'MaxSpeedup',max(s,[],'omitnan'), ...
                    'MedianSeconds',median(R.MedianOptimizedSeconds,'omitnan')); %#ok<AGROW>
            end
        end
    end
    if isempty(rows), A=table(); else, A=struct2table(vertcat(rows{:})); end
end


function A=aggregate_auto_by_M(T)
    Ms=unique(T.M(:).'); rows={}; q=0;
    for M=Ms
        R=T(T.M==M,:); s=R.ReferenceOverOptimizedSpeedup;
        q=q+1; rows{q,1}=struct('M',M,'NScenarios',height(R), ...
            'MedianSpeedup',median(s,'omitnan'),'MeanSpeedup',mean(s,'omitnan'), ...
            'MinSpeedup',min(s,[],'omitnan'),'MaxSpeedup',max(s,[],'omitnan')); %#ok<AGROW>
    end
    if isempty(rows), A=table(); else, A=struct2table(vertcat(rows{:})); end
end


function [recommended,bestMedian]=recommend_by_M(A,threshold)
    Ms=unique(A.M(:).'); rec={}; bm={}; qr=0; qb=0;
    for M=Ms
        R=A(A.M==M,:);
        % Highest median, irrespective of worst-case regression.
        maxMed=max(R.MedianSpeedup); cand=R(abs(R.MedianSpeedup-maxMed)<1e-12,:);
        [~,ii]=max(cand.MinSpeedup); B=cand(ii,:);
        qb=qb+1; bm{qb,1}=table_row_to_struct(B,false,threshold); %#ok<AGROW>

        eligible=R(R.MinSpeedup>=threshold,:);
        if ~isempty(eligible)
            maxEligibleMedian=max(eligible.MedianSpeedup);
            cand=eligible(abs(eligible.MedianSpeedup-maxEligibleMedian)<1e-12,:);
            [~,ii]=max(cand.MinSpeedup); S=cand(ii,:); noRegression=true;
        else
            maxMin=max(R.MinSpeedup); cand=R(abs(R.MinSpeedup-maxMin)<1e-12,:);
            [~,ii]=max(cand.MedianSpeedup); S=cand(ii,:); noRegression=false;
        end
        qr=qr+1; rec{qr,1}=table_row_to_struct(S,noRegression,threshold); %#ok<AGROW>
    end
    recommended=struct2table(vertcat(rec{:})); bestMedian=struct2table(vertcat(bm{:}));
end


function s=table_row_to_struct(R,noRegression,threshold)
    s=struct('M',R.M(1),'BlockSize',R.BlockSize(1),'BatchSize',R.BatchSize(1), ...
        'NScenarios',R.NScenarios(1),'MedianSpeedup',R.MedianSpeedup(1), ...
        'MeanSpeedup',R.MeanSpeedup(1),'MinSpeedup',R.MinSpeedup(1), ...
        'MaxSpeedup',R.MaxSpeedup(1),'MedianSeconds',R.MedianSeconds(1), ...
        'MeetsNoRegressionThreshold',logical(noRegression), ...
        'NoRegressionThreshold',threshold);
end


function warm_eval(x,cfg,dt,scale,block,nWarmup)
    for k=1:nWarmup %#ok<NASGU>
        AMI_noCSI_fast_logh_adaptive_grid(x,cfg.px,cfg,dt,8,scale,block);
    end
end
function [v,tmed,tmean,tmin]=timed_eval(x,cfg,dt,scale,block,nRepeats)
    times=nan(nRepeats,1); vals=nan(nRepeats,1);
    for r=1:nRepeats
        tt=tic; vals(r)=AMI_noCSI_fast_logh_adaptive_grid(x,cfg.px,cfg,dt,8,scale,block); times(r)=toc(tt);
    end
    v=vals(end); tmed=median(times); tmean=mean(times); tmin=min(times);
end
function assert_agreement(v,ref,tol,M,snr,sig,geometry,label)
    err=abs(v-ref);
    if err>tol*max(1,abs(ref))
        error('diagnose_i10p2p1_runtime_tuning_sweep:NumericalMismatch', ...
            'M=%d SNR=%.1f sigma=%.3f geometry=%s %s mismatch %.3e.', ...
            M,snr,sig,char(geometry),label,err);
    end
end
function x=make_geometry(mode,M,Pavg,dmin)
    switch string(mode)
        case "uniform"
            x=define_constellation(M,Pavg,true,"mean");
        case "wide"
            x=wide_geometry(M,Pavg,dmin);
        otherwise
            error('diagnose_i10p2p1_runtime_tuning_sweep:BadGeometry','Unknown geometry %s.',mode);
    end
end
function x=wide_geometry(M,Pavg,dmin)
    K=max(1,ceil(M/8)); b=(0:M-1).'*dmin; surplus=Pavg-mean(b);
    if surplus<0,error('diagnose_i10p2p1_runtime_tuning_sweep:Infeasible','Infeasible wide geometry.');end
    q=zeros(M,1); q(end-K+1:end)=surplus*M/K; x=b+q;
end
function tf=cases_valid(v)
    tf=isnumeric(v)&&ismatrix(v)&&~isempty(v)&&any(size(v,2)==[3 4])&&all(isfinite(v(:)))&& ...
        all(v(:,1)>=2)&&all(mod(v(:,1),1)==0)&&all(v(:,3)>=0);
    if tf && size(v,2)==4, tf=all(v(:,4)>0); end
end
function tf=geometry_modes_valid(v)
    try
        s=lower(string(v)); tf=~isempty(s)&&all(ismember(s,["uniform","wide"]));
    catch
        tf=false;
    end
end
function tf=positive_integer_vector(v),tf=isnumeric(v)&&isvector(v)&&~isempty(v)&&all(isfinite(v))&&all(v>=1)&&all(mod(v,1)==0);end
function tf=positive_integer(v),tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>=1&&mod(v,1)==0;end
function tf=nonnegative_integer(v),tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>=0&&mod(v,1)==0;end
function tf=positive_scalar(v),tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&isreal(v)&&v>0;end
function tf=nonnegative_scalar(v),tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&isreal(v)&&v>=0;end
function t=sanitize_label(v),t=strtrim(char(string(v)));if isempty(t),t='run';return;end,t=regexprep(t,'[^A-Za-z0-9._-]+','_');t=regexprep(t,'_+','_');t=regexprep(t,'^[._-]+|[._-]+$','');if isempty(t),t='run';end,end
function [d,id]=make_unique_dir(parent,id),if ~exist(parent,'dir'),mkdir(parent);end,d=fullfile(parent,id);base=id;n=1;while exist(d,'dir'),n=n+1;id=sprintf('%s_%02d',base,n);d=fullfile(parent,id);end,[ok,msg]=mkdir(d);if ~ok,error('diagnose_i10p2p1_runtime_tuning_sweep:mkdir','%s',msg);end,end
