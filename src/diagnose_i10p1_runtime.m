function d = diagnose_i10p1_runtime(varargin)
%DIAGNOSE_I10P1_RUNTIME  Benchmark the I.10.1 implementation-only speedups.
%
% Representative defaults:
%   (4,20,.1)   easy/base-dt case
%   (16,30,.2)  I.9 boundary case
%   (32,30,.3)  heavy refined-dt case
%
% For each case a wide feasible constellation is used so the adaptive y-grid
% resembles the numerically expensive stress geometries. The diagnostic:
%   1) times a private copy of the pre-I.10.1 legacy evaluator at block=512;
%   2) times the optimized cached evaluator at several block sizes;
%   3) verifies numerical agreement against the legacy value;
%   4) reports per-case speedup and the fastest tested block size.
%
% This script does not run SA and does not change production settings.

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'Cases',[4 20 0.1;16 30 0.2;32 30 0.3],@cases_valid);
    addParameter(p,'BlockSizes',[256 512 1024 2048],@positive_integer_vector);
    addParameter(p,'nWarmup',1,@nonnegative_integer);
    addParameter(p,'nRepeats',3,@positive_integer);
    addParameter(p,'AdaptiveScale',1.1,@(v) isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>=1);
    addParameter(p,'minGap',0.01,@nonnegative_scalar);
    addParameter(p,'AbsAgreementTolerance',5e-12,@positive_scalar);
    addParameter(p,'WriteCSV',true,@(v) islogical(v)&&isscalar(v));
    addParameter(p,'RunLabel','',@(v) ischar(v)||(isstring(v)&&isscalar(v)));
    addParameter(p,'ResultsRoot',fso_result_utils.default_results_root(),@(v) ischar(v)||isstring(v));
    parse(p,varargin{:}); o=p.Results;

    C=double(o.Cases); blocks=unique(double(o.BlockSizes(:).'),'stable');
    parent=fullfile(char(o.ResultsRoot),'i10p1_runtime');
    label=sanitize_label(o.RunLabel);
    utc=char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmssSSS''Z'''));
    [outDir,runId]=make_unique_dir(parent,sprintf('%s_%s',label,utc));

    fprintf('\n============================================================\n');
    fprintf('Commit I.10.1 - evaluator runtime benchmark\n');
    fprintf('Cases [M SNR sigma_X^2]:\n'); disp(C);
    fprintf('Block sizes: %s | repeats=%d | warmup=%d\n',mat2str(blocks),o.nRepeats,o.nWarmup);
    fprintf('Results: %s\n',outDir);
    fprintf('============================================================\n\n');

    rows={}; q=0;
    for k=1:size(C,1)
        M=C(k,1); snr=C(k,2); sig=C(k,3); [dt,meta]=select_logh_dt(M,snr,sig);
        cfg=build_fso_config(M,1,snr,sig,'minGap',o.minGap,'pinZero',true, ...
            'saMaxIter',1,'saNStarts',1,'YGridMode','candidate-adaptive', ...
            'YGridScale',o.AdaptiveScale,'loghDt',dt,'loghSpanSigma',8,'loghYBlockSize',512);
        x=wide_geometry(M,1,o.minGap); AMI_functions.assert_constellation_feasible(x,cfg,1e-10);

        t=tic; legacy=legacy_adaptive_eval(x,cfg.px,cfg,dt,8,o.AdaptiveScale,512); legacySec=toc(t);

        for b=blocks
            for w=1:o.nWarmup %#ok<NASGU>
                AMI_noCSI_fast_logh_adaptive_grid(x,cfg.px,cfg,dt,8,o.AdaptiveScale,b);
            end
            times=nan(o.nRepeats,1); vals=nan(o.nRepeats,1);
            for r=1:o.nRepeats
                tt=tic;
                vals(r)=AMI_noCSI_fast_logh_adaptive_grid(x,cfg.px,cfg,dt,8,o.AdaptiveScale,b);
                times(r)=toc(tt);
            end
            v=vals(end); err=abs(v-legacy);
            if err>o.AbsAgreementTolerance*max(1,abs(legacy))
                error('diagnose_i10p1_runtime:NumericalMismatch', ...
                    'M=%d SNR=%.1f sigma=%.3f block=%d mismatch %.3e.',M,snr,sig,b,err);
            end
            q=q+1;
            rows(q,:)={M,sig,snr,dt,logical(meta.isRefined),max(x),b,legacy,legacySec, ...
                v,median(times),mean(times),min(times),legacySec/median(times),err}; %#ok<AGROW>
        end
        fprintf('  completed %d/%d physical cases\n',k,size(C,1));
    end

    T=cell2table(rows,'VariableNames',{'M','SigmaX2','SNRdB','Dt','IsRefined','XMax','BlockSize', ...
        'LegacyAMI','LegacySeconds','OptimizedAMI','MedianOptimizedSeconds','MeanOptimizedSeconds', ...
        'MinOptimizedSeconds','LegacyOverOptimizedSpeedup','AbsDifference'});
    numNames=setdiff(T.Properties.VariableNames,{'IsRefined'});
    for j=1:numel(numNames), if iscell(T.(numNames{j})),T.(numNames{j})=cell2mat(T.(numNames{j}));end,end
    if iscell(T.IsRefined),T.IsRefined=logical(cell2mat(T.IsRefined));end

    bestRows=false(height(T),1);
    for k=1:size(C,1)
        mask=T.M==C(k,1)&T.SNRdB==C(k,2)&abs(T.SigmaX2-C(k,3))<1e-12;
        idx=find(mask); [~,ii]=min(T.MedianOptimizedSeconds(idx)); bestRows(idx(ii))=true;
    end
    T.IsFastestTestedBlock=bestRows;

    B=T(T.BlockSize==512,:);
    summary=struct();
    summary.meanSpeedupAt512=mean(B.LegacyOverOptimizedSpeedup,'omitnan');
    summary.medianSpeedupAt512=median(B.LegacyOverOptimizedSpeedup,'omitnan');
    summary.minSpeedupAt512=min(B.LegacyOverOptimizedSpeedup,[],'omitnan');
    summary.maxSpeedupAt512=max(B.LegacyOverOptimizedSpeedup,[],'omitnan');
    summary.maxAbsDifference=max(T.AbsDifference,[],'omitnan');
    summary.bestByCase=T(T.IsFastestTestedBlock,:);

    d=struct(); d.run=struct('id',runId,'outputDirectory',outDir); d.results=T; d.summary=summary;
    save(fullfile(outDir,'i10p1Runtime.mat'),'d','-v7.3');
    if o.WriteCSV,writetable(T,fullfile(outDir,'i10p1Runtime.csv'));end

    fprintf('\nI.10.1 RUNTIME SUMMARY\n');
    fprintf('block=512 legacy/optimized speedup mean/median : %.3fx / %.3fx\n', ...
        summary.meanSpeedupAt512,summary.medianSpeedupAt512);
    fprintf('block=512 speedup min/max                     : %.3fx / %.3fx\n', ...
        summary.minSpeedupAt512,summary.maxSpeedupAt512);
    fprintf('max numerical |optimized-legacy|              : %.3e\n',summary.maxAbsDifference);
    fprintf('\nFastest tested block by case:\n');
    disp(summary.bestByCase(:,{'M','SigmaX2','SNRdB','Dt','XMax','BlockSize','MedianOptimizedSeconds','LegacyOverOptimizedSpeedup'}));
end


function mi=legacy_adaptive_eval(x,px,params,dt,spanSigma,scale,block)
    y=AMI_functions.build_noCSI_y_grid(params,scale*max(abs(x)));
    mi=legacy_logh_grid(x,px,params,dt,spanSigma,y,block);
end
function mi_bits=legacy_logh_grid(x,px,params,dt,spanSigma,y_grid,yBlockSize)
    x=real(x(:).'); px=real(px(:).'); px=px/sum(px); y_grid=real(y_grid(:).');
    M=numel(x); Ny=numel(y_grid);
    if params.sig_t<1e-12
        means=params.R*x(:); d2=bsxfun(@minus,y_grid,means).^2;
        logPyx=-0.5*log(2*pi*params.sigma_n_sq)-d2/(2*params.sigma_n_sq);
        mi_bits=legacy_mi(logPyx,px,y_grid); return;
    end
    tMin=params.mu_t-spanSigma*params.sig_t; tMax=params.mu_t+spanSigma*params.sig_t;
    tGrid=tMin:dt:tMax;
    if isempty(tGrid)||tGrid(end)<tMax-100*eps(max(1,abs(tMax))),tGrid=[tGrid tMax]; %#ok<AGROW>
    else,tGrid(end)=tMax;end
    tGrid=unique(tGrid,'stable'); qWeights=legacy_trap_weights(tGrid);
    fT=exp(-0.5*((tGrid-params.mu_t)/params.sig_t).^2)/(sqrt(2*pi)*params.sig_t);
    q=qWeights.*fT; keep=q>0&isfinite(q); tGrid=tGrid(keep); q=q(keep);
    logQ=log(q(:)); hVals=exp(tGrid(:)); inv2s2=1/(2*params.sigma_n_sq);
    logNorm=-0.5*log(2*pi*params.sigma_n_sq); logPyx=zeros(M,Ny);
    for j=1:M
        means=params.R*hVals*x(j);
        for b0=1:yBlockSize:Ny
            b1=min(Ny,b0+yBlockSize-1); yb=y_grid(b0:b1);
            d2=bsxfun(@minus,yb,means).^2;
            logTerms=bsxfun(@plus,logQ+logNorm,-inv2s2*d2);
            logPyx(j,b0:b1)=legacy_lse(logTerms,1);
        end
    end
    mi_bits=legacy_mi(logPyx,px,y_grid);
end
function mi=legacy_mi(logPyx,px,y)
    logPxPyx=bsxfun(@plus,log(px(:)),logPyx); logPy=legacy_lse(logPxPyx,1);
    M=size(logPyx,1); MI=0;
    for j=1:M
        Pyxj=exp(logPyx(j,:)); lr=(logPyx(j,:)-logPy)/log(2);
        z=Pyxj.*lr; z(~isfinite(z))=0; MI=MI+px(j)*trapz(y,z);
    end
    mi=max(0,MI);
end
function w=legacy_trap_weights(x)
    x=x(:).'; n=numel(x); if n==1,w=1;return;end
    dx=diff(x); w=zeros(1,n); w(1)=dx(1)/2; w(end)=dx(end)/2;
    if n>2,w(2:end-1)=(dx(1:end-1)+dx(2:end))/2;end
end
function s=legacy_lse(A,dim),m=max(A,[],dim);s=m+log(sum(exp(bsxfun(@minus,A,m)),dim));end

function x=wide_geometry(M,Pavg,dmin)
    K=max(1,ceil(M/8)); b=(0:M-1).'*dmin; surplus=Pavg-mean(b);
    if surplus<0,error('diagnose_i10p1_runtime:Infeasible','Infeasible wide geometry.');end
    q=zeros(M,1); q(end-K+1:end)=surplus*M/K; x=b+q;
end
function tf=cases_valid(v),tf=isnumeric(v)&&ismatrix(v)&&size(v,2)==3&&~isempty(v)&&all(isfinite(v(:)))&&all(v(:,1)>=2)&&all(mod(v(:,1),1)==0)&&all(v(:,3)>=0);end
function tf=positive_integer_vector(v),tf=isnumeric(v)&&isvector(v)&&~isempty(v)&&all(isfinite(v))&&all(v>=1)&&all(mod(v,1)==0);end
function tf=positive_integer(v),tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>=1&&mod(v,1)==0;end
function tf=nonnegative_integer(v),tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>=0&&mod(v,1)==0;end
function tf=positive_scalar(v),tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&isreal(v)&&v>0;end
function tf=nonnegative_scalar(v),tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&isreal(v)&&v>=0;end
function t=sanitize_label(v),t=strtrim(char(string(v)));if isempty(t),t='run';return;end,t=regexprep(t,'[^A-Za-z0-9._-]+','_');t=regexprep(t,'_+','_');t=regexprep(t,'^[._-]+|[._-]+$','');if isempty(t),t='run';end,end
function [d,id]=make_unique_dir(parent,id),if ~exist(parent,'dir'),mkdir(parent);end,d=fullfile(parent,id);base=id;n=1;while exist(d,'dir'),n=n+1;id=sprintf('%s_%02d',base,n);d=fullfile(parent,id);end,[ok,msg]=mkdir(d);if ~ok,error('diagnose_i10p1_runtime:mkdir','%s',msg);end,end
