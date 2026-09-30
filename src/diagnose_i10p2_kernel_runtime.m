function d = diagnose_i10p2_kernel_runtime(varargin)
%DIAGNOSE_I10P2_KERNEL_RUNTIME  Benchmark Commit-I.10.2 kernel choices.
%
% Compares the exact I.10.1 cached kernel (fixed block=512, one symbol at a
% time, historical unique-grid path) against:
%   1) the default I.10.2 auto runtime policy;
%   2) a block-size x symbol-batch sweep using the I.10.2 zero-level fast path.
%
% No SA is run and no scientific settings are changed.

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'Cases',[4 20 0.1;8 25 0.2;16 30 0.2;32 30 0.3],@cases_valid);
    addParameter(p,'BlockSizes',[128 256 512 1024 2048],@positive_integer_vector);
    addParameter(p,'BatchSizes',[1 2 4],@positive_integer_vector);
    addParameter(p,'nWarmup',1,@nonnegative_integer);
    addParameter(p,'nRepeats',3,@positive_integer);
    addParameter(p,'AdaptiveScale',1.1,@(v) isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>=1);
    addParameter(p,'minGap',0.01,@nonnegative_scalar);
    addParameter(p,'AbsAgreementTolerance',5e-12,@positive_scalar);
    addParameter(p,'WriteCSV',true,@(v) islogical(v)&&isscalar(v));
    addParameter(p,'RunLabel','',@(v) ischar(v)||(isstring(v)&&isscalar(v)));
    addParameter(p,'ResultsRoot',fso_result_utils.default_results_root(),@(v) ischar(v)||isstring(v));
    parse(p,varargin{:}); o=p.Results;

    C=double(o.Cases); blocks=unique(double(o.BlockSizes(:).'),'stable'); batches=unique(double(o.BatchSizes(:).'),'stable');
    parent=fullfile(char(o.ResultsRoot),'i10p2_kernel_runtime');
    label=sanitize_label(o.RunLabel); utc=char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmssSSS''Z'''));
    [outDir,runId]=make_unique_dir(parent,sprintf('%s_%s',label,utc));

    fprintf('\n============================================================\n');
    fprintf('Commit I.10.2 - likelihood-kernel runtime benchmark\n');
    fprintf('Cases [M SNR sigma_X^2]:\n'); disp(C);
    fprintf('Blocks=%s | batches=%s | repeats=%d | warmup=%d\n', ...
        mat2str(blocks),mat2str(batches),o.nRepeats,o.nWarmup);
    fprintf('Results: %s\n',outDir);
    fprintf('============================================================\n\n');

    rows={}; q=0;
    for k=1:size(C,1)
        M=C(k,1); snr=C(k,2); sig=C(k,3); [dt,metaDt]=select_logh_dt(M,snr,sig);
        cfg=build_fso_config(M,1,snr,sig,'minGap',o.minGap,'pinZero',true, ...
            'saMaxIter',1,'saNStarts',1,'YGridMode','candidate-adaptive', ...
            'YGridScale',o.AdaptiveScale,'loghDt',dt,'loghSpanSigma',8,'loghYBlockSize',512);
        x=wide_geometry(M,1,o.minGap); AMI_functions.assert_constellation_feasible(x,cfg,1e-10);

        ref=cfg; ref.loghRuntimeKernel='i10p1'; ref.loghRuntimeBlockPolicy='fixed';
        ref.loghSymbolBatchSize=1; ref.loghZeroFastPath=false;
        for w=1:o.nWarmup %#ok<NASGU>
            AMI_noCSI_fast_logh_adaptive_grid(x,ref.px,ref,dt,8,o.AdaptiveScale,512);
        end
        refTimes=nan(o.nRepeats,1); refVals=nan(o.nRepeats,1);
        for r=1:o.nRepeats
            tt=tic; refVals(r)=AMI_noCSI_fast_logh_adaptive_grid(x,ref.px,ref,dt,8,o.AdaptiveScale,512); refTimes(r)=toc(tt);
        end
        refAMI=refVals(end); refSec=median(refTimes);

        % Default I.10.2 auto policy.
        auto=cfg;
        [~,dAuto]=AMI_noCSI_fast_logh_adaptive_grid(x,auto.px,auto,dt,8,o.AdaptiveScale,512);
        for w=1:o.nWarmup %#ok<NASGU>
            AMI_noCSI_fast_logh_adaptive_grid(x,auto.px,auto,dt,8,o.AdaptiveScale,512);
        end
        [v,tmed,tmean,tmin]=timed_eval(x,auto,dt,o.AdaptiveScale,512,o.nRepeats);
        assert_agreement(v,refAMI,o.AbsAgreementTolerance,M,snr,sig,'auto');
        q=q+1; rows(q,:)={"auto",M,sig,snr,dt,logical(metaDt.isRefined),max(x),512, ...
            dAuto.yBlockSizeEffective,dAuto.symbolBatchSize,refAMI,refSec,v,tmed,tmean,tmin,refSec/tmed,abs(v-refAMI)}; %#ok<AGROW>

        % Explicit implementation-only sweep.
        for b=blocks
            for sb=batches
                tuned=cfg; tuned.loghRuntimeKernel='i10p2'; tuned.loghRuntimeBlockPolicy='fixed';
                tuned.loghSymbolBatchSize=sb; tuned.loghZeroFastPath=true;
                [~,dTune]=AMI_noCSI_fast_logh_adaptive_grid(x,tuned.px,tuned,dt,8,o.AdaptiveScale,b);
                for w=1:o.nWarmup %#ok<NASGU>
                    AMI_noCSI_fast_logh_adaptive_grid(x,tuned.px,tuned,dt,8,o.AdaptiveScale,b);
                end
                [v,tmed,tmean,tmin]=timed_eval(x,tuned,dt,o.AdaptiveScale,b,o.nRepeats);
                assert_agreement(v,refAMI,o.AbsAgreementTolerance,M,snr,sig,sprintf('b%d-s%d',b,sb));
                q=q+1; rows(q,:)={"tuned",M,sig,snr,dt,logical(metaDt.isRefined),max(x),b, ...
                    dTune.yBlockSizeEffective,dTune.symbolBatchSize,refAMI,refSec,v,tmed,tmean,tmin,refSec/tmed,abs(v-refAMI)}; %#ok<AGROW>
            end
        end
        fprintf('  completed %d/%d physical cases\n',k,size(C,1));
    end

    T=cell2table(rows,'VariableNames',{'Mode','M','SigmaX2','SNRdB','Dt','IsRefined','XMax', ...
        'BlockRequested','BlockEffective','BatchSize','ReferenceAMI','ReferenceSeconds','OptimizedAMI', ...
        'MedianOptimizedSeconds','MeanOptimizedSeconds','MinOptimizedSeconds','ReferenceOverOptimizedSpeedup','AbsDifference'});
    T.Mode=string(T.Mode);
    if iscell(T.IsRefined)
        T.IsRefined=logical(cell2mat(T.IsRefined));
    else
        T.IsRefined=logical(T.IsRefined);
    end
    numericNames=setdiff(T.Properties.VariableNames,{'Mode','IsRefined'});
    for j=1:numel(numericNames), if iscell(T.(numericNames{j})),T.(numericNames{j})=cell2mat(T.(numericNames{j}));end,end

    best=false(height(T),1);
    for k=1:size(C,1)
        mask=T.M==C(k,1)&T.SNRdB==C(k,2)&abs(T.SigmaX2-C(k,3))<1e-12;
        idx=find(mask); [~,ii]=min(T.MedianOptimizedSeconds(idx)); best(idx(ii))=true;
    end
    T.IsFastestTested=best;

    autoRows=T(T.Mode=="auto",:); bestRows=T(T.IsFastestTested,:);
    summary=struct();
    summary.autoByCase=autoRows;
    summary.bestByCase=bestRows;
    summary.meanAutoSpeedup=mean(autoRows.ReferenceOverOptimizedSpeedup,'omitnan');
    summary.medianAutoSpeedup=median(autoRows.ReferenceOverOptimizedSpeedup,'omitnan');
    summary.minAutoSpeedup=min(autoRows.ReferenceOverOptimizedSpeedup,[],'omitnan');
    summary.maxAutoSpeedup=max(autoRows.ReferenceOverOptimizedSpeedup,[],'omitnan');
    summary.maxAbsDifference=max(T.AbsDifference,[],'omitnan');

    d=struct(); d.run=struct('id',runId,'outputDirectory',outDir); d.results=T; d.summary=summary;
    save(fullfile(outDir,'i10p2KernelRuntime.mat'),'d','-v7.3');
    if o.WriteCSV,writetable(T,fullfile(outDir,'i10p2KernelRuntime.csv'));end

    fprintf('\nI.10.2 KERNEL RUNTIME SUMMARY\n');
    fprintf('auto speedup mean/median/min/max : %.3fx / %.3fx / %.3fx / %.3fx\n', ...
        summary.meanAutoSpeedup,summary.medianAutoSpeedup,summary.minAutoSpeedup,summary.maxAutoSpeedup);
    fprintf('max numerical |optimized-I10.1| : %.3e\n',summary.maxAbsDifference);
    fprintf('\nAuto policy by case:\n');
    disp(autoRows(:,{'M','SigmaX2','SNRdB','Dt','BlockEffective','BatchSize','MedianOptimizedSeconds','ReferenceOverOptimizedSpeedup'}));
    fprintf('\nFastest tested combination by case:\n');
    disp(bestRows(:,{'M','SigmaX2','SNRdB','Dt','BlockEffective','BatchSize','MedianOptimizedSeconds','ReferenceOverOptimizedSpeedup'}));
end


function [v,tmed,tmean,tmin]=timed_eval(x,cfg,dt,scale,block,nRepeats)
    times=nan(nRepeats,1); vals=nan(nRepeats,1);
    for r=1:nRepeats
        tt=tic; vals(r)=AMI_noCSI_fast_logh_adaptive_grid(x,cfg.px,cfg,dt,8,scale,block); times(r)=toc(tt);
    end
    v=vals(end); tmed=median(times); tmean=mean(times); tmin=min(times);
end
function assert_agreement(v,ref,tol,M,snr,sig,label)
    err=abs(v-ref);
    if err>tol*max(1,abs(ref))
        error('diagnose_i10p2_kernel_runtime:NumericalMismatch', ...
            'M=%d SNR=%.1f sigma=%.3f %s mismatch %.3e.',M,snr,sig,label,err);
    end
end
function x=wide_geometry(M,Pavg,dmin)
    K=max(1,ceil(M/8)); b=(0:M-1).'*dmin; surplus=Pavg-mean(b);
    if surplus<0,error('diagnose_i10p2_kernel_runtime:Infeasible','Infeasible wide geometry.');end
    q=zeros(M,1); q(end-K+1:end)=surplus*M/K; x=b+q;
end
function tf=cases_valid(v),tf=isnumeric(v)&&ismatrix(v)&&size(v,2)==3&&~isempty(v)&&all(isfinite(v(:)))&&all(v(:,1)>=2)&&all(mod(v(:,1),1)==0)&&all(v(:,3)>=0);end
function tf=positive_integer_vector(v),tf=isnumeric(v)&&isvector(v)&&~isempty(v)&&all(isfinite(v))&&all(v>=1)&&all(mod(v,1)==0);end
function tf=positive_integer(v),tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>=1&&mod(v,1)==0;end
function tf=nonnegative_integer(v),tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>=0&&mod(v,1)==0;end
function tf=positive_scalar(v),tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&isreal(v)&&v>0;end
function tf=nonnegative_scalar(v),tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&isreal(v)&&v>=0;end
function t=sanitize_label(v),t=strtrim(char(string(v)));if isempty(t),t='run';return;end,t=regexprep(t,'[^A-Za-z0-9._-]+','_');t=regexprep(t,'_+','_');t=regexprep(t,'^[._-]+|[._-]+$','');if isempty(t),t='run';end,end
function [d,id]=make_unique_dir(parent,id),if ~exist(parent,'dir'),mkdir(parent);end,d=fullfile(parent,id);base=id;n=1;while exist(d,'dir'),n=n+1;id=sprintf('%s_%02d',base,n);d=fullfile(parent,id);end,[ok,msg]=mkdir(d);if ~ok,error('diagnose_i10p2_kernel_runtime:mkdir','%s',msg);end,end
