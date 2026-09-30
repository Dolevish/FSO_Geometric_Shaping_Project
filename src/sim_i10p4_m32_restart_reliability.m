function d = sim_i10p4_m32_restart_reliability(varargin)
%SIM_I10P4_M32_RESTART_RELIABILITY  Focused 4-vs-6 restart SA pilot.
%
% Commit I.10.4 asks whether the production candidate budget of four SA
% restarts is sufficient for the two M=32 high-SNR cases that showed the
% largest between-replicate variability in I.10.  The experiment runs ONE
% six-restart optimization per replicate, then compares the validated winner
% among starts 1:4 with the validated winner among starts 1:6.  Because the
% first four starts are literally a prefix of the six-start run, this is a
% paired/nested comparison with identical initialization and RNG substreams.
%
% Default scientific settings:
%   M=32, SNR=30 dB, sigma_X^2=[0.2 0.3]
%   dt selected by the frozen I.10.3 case-adaptive policy (=0.005 here)
%   candidate budget = 4 restarts, extended budget = 6 restarts
%   4000 iterations/restart, 3 independent replicates
%   T0=0.4, Tf=1e-3, BaseStd0=0.15, ItersPerTemp=50
%
% Hard decision criterion for restart count:
%   max(best_6_validated - best_4_validated) <= ExtraRestartTolerance
% with default tolerance 1e-3 bit/symbol.  Repeatability across replicates is
% reported separately and is NOT assigned an arbitrary hard threshold.

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'M',32,@positive_integer);
    addParameter(p,'SNRdB',30,@finite_scalar);
    addParameter(p,'TurbulenceVec',[0.2 0.3],@nonnegative_vector);
    addParameter(p,'P_avg',1,@positive_scalar);
    addParameter(p,'minGap',0.01,@nonnegative_scalar);
    addParameter(p,'pinZero',true,@(v)islogical(v)&&isscalar(v));
    addParameter(p,'CandidateStarts',4,@(v)positive_integer(v)&&v>=2);
    addParameter(p,'ExtendedStarts',6,@(v)positive_integer(v)&&v>=3);
    addParameter(p,'nReplicates',3,@positive_integer);
    addParameter(p,'saMaxIter',4000,@positive_integer);
    addParameter(p,'T0',0.4,@positive_scalar);
    addParameter(p,'Tf',1e-3,@positive_scalar);
    addParameter(p,'BaseStd0',0.15,@positive_scalar);
    addParameter(p,'ItersPerTemp',50,@positive_integer);
    addParameter(p,'baseSeed',20260928,@nonnegative_scalar);
    addParameter(p,'historyEvery',50,@positive_integer);
    addParameter(p,'ExtraRestartTolerance',1e-3,@positive_scalar);
    addParameter(p,'FastValidationTolerance',1e-3,@positive_scalar);
    addParameter(p,'UseParallelCases',true,@(v)islogical(v)&&isscalar(v));
    addParameter(p,'WriteCSV',true,@(v)islogical(v)&&isscalar(v));
    addParameter(p,'RunLabel','',@(v)ischar(v)||(isstring(v)&&isscalar(v)));
    addParameter(p,'ResultsRoot',fso_result_utils.default_results_root(),@(v)ischar(v)||isstring(v));
    parse(p,varargin{:}); o=p.Results;

    if o.CandidateStarts>=o.ExtendedStarts
        error('sim_i10p4_m32_restart_reliability:RestartOrder', ...
            'CandidateStarts must be smaller than ExtendedStarts.');
    end
    if ~(o.T0>o.Tf)
        error('sim_i10p4_m32_restart_reliability:TemperatureOrder','T0 must exceed Tf.');
    end

    sigVec=unique(double(o.TurbulenceVec(:).'),'stable');
    parent=fullfile(char(o.ResultsRoot),'i10p4_m32_restart_reliability');
    label=sanitize_label(o.RunLabel);
    utc=char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmssSSS''Z'''));
    token=sprintf('%s_M%d_SNR%s_I%d_C%d_E%d_R%d_%s',label,o.M,num_token(o.SNRdB,2), ...
        o.saMaxIter,o.CandidateStarts,o.ExtendedStarts,o.nReplicates,utc);
    [outDir,runId]=make_unique_dir(parent,token);
    prodRoot=fullfile(outDir,'production'); mkdir(prodRoot);

    fprintf('\n============================================================\n');
    fprintf('Commit I.10.4 - focused restart reliability pilot\n');
    fprintf('M=%d | SNR=%.1f dB | sigma_X^2=%s\n',o.M,o.SNRdB,mat2str(sigVec));
    fprintf('Nested restart comparison: first %d vs all %d starts\n',o.CandidateStarts,o.ExtendedStarts);
    fprintf('SA=%d iterations/restart | replicates=%d | extra-start tolerance=%.3e bit/symbol\n', ...
        o.saMaxIter,o.nReplicates,o.ExtraRestartTolerance);
    fprintf('Results: %s\n',outDir);
    fprintf('============================================================\n\n');

    b=sim_revision_AMI_vs_SNR( ...
        'MVec',o.M,'SNRVec',o.SNRdB,'TurbulenceVec',sigVec,'P_avg',o.P_avg, ...
        'minGap',o.minGap,'pinZero',o.pinZero,'FastFadingMethod','logh', ...
        'LoghDtPolicy','case-adaptive','loghDt',0.01,'loghDtRefined',0.005, ...
        'loghSpanSigma',8,'loghYBlockSize',512,'YGridMode','candidate-adaptive','YGridScale',1.1, ...
        'saMaxIter',o.saMaxIter,'saNStarts',o.ExtendedStarts,'nReplicates',o.nReplicates, ...
        'T0',o.T0,'Tf',o.Tf,'BaseStd0',o.BaseStd0,'ItersPerTemp',o.ItersPerTemp, ...
        'baseSeed',o.baseSeed,'historyEvery',o.historyEvery, ...
        'UseParallelCases',o.UseParallelCases,'MakePlots',false,'SaveFigures',false, ...
        'ResultsRoot',prodRoot,'RunLabel','I10p4_nested6');

    nSig=numel(sigVec); nRep=double(o.nReplicates); nRows=nSig*nRep;
    SigmaX2=nan(nRows,1); Replicate=nan(nRows,1); Seed=nan(nRows,1);
    SelectedDt=nan(nRows,1); IsRefined=false(nRows,1); RuleVersion=strings(nRows,1);
    BestCandidateAMI=nan(nRows,1); BestExtendedAMI=nan(nRows,1); ExtraValidatedGain=nan(nRows,1);
    BestCandidateStart=nan(nRows,1); BestExtendedStart=nan(nRows,1);
    ExtraStartWins=false(nRows,1); MaterialExtraGain=false(nRows,1);
    CandidateK99=nan(nRows,1); ExtendedK99=nan(nRows,1);
    CandidateRestartStd=nan(nRows,1); ExtendedRestartStd=nan(nRows,1);
    CandidateMaxFVGap=nan(nRows,1); ExtendedMaxFVGap=nan(nRows,1);
    CandidateSelectionChanged=false(nRows,1); ExtendedSelectionChanged=false(nRows,1);
    CandidateSerialSeconds=nan(nRows,1); ExtendedSerialSeconds=nan(nRows,1);
    ExtendedOverCandidateRuntimeRatio=nan(nRows,1); MeasuredExtendedWallSeconds=nan(nRows,1);
    MaxAbsConstellationDelta=nan(nRows,1); CaseFile=strings(nRows,1);

    q=0;
    for iSig=1:nSig
        for r=1:nRep
            q=q+1; path=b.caseFiles{1,iSig,1,r}; s=load(path,'c'); c=s.c; starts=c.starts;
            if numel(starts)~=o.ExtendedStarts
                error('sim_i10p4_m32_restart_reliability:StartCount', ...
                    'Expected %d starts, found %d in %s.',o.ExtendedStarts,numel(starts),path);
            end
            vals=[starts.bestMIValidated]; fast=[starts.bestMIFast];
            [vC,iC]=max(vals(1:o.CandidateStarts)); [vE,iE]=max(vals(1:o.ExtendedStarts));
            [~,iFastC]=max(fast(1:o.CandidateStarts)); [~,iFastE]=max(fast(1:o.ExtendedStarts));

            SigmaX2(q)=sigVec(iSig); Replicate(q)=r; Seed(q)=c.case.seed;
            SelectedDt(q)=c.integration.loghDtSelected; IsRefined(q)=c.integration.loghDtIsRefined;
            RuleVersion(q)=string(c.integration.loghDtRuleVersion);
            BestCandidateAMI(q)=vC; BestExtendedAMI(q)=vE; ExtraValidatedGain(q)=vE-vC;
            BestCandidateStart(q)=iC; BestExtendedStart(q)=iE;
            ExtraStartWins(q)=iE>o.CandidateStarts;
            MaterialExtraGain(q)=ExtraValidatedGain(q)>o.ExtraRestartTolerance;
            CandidateK99(q)=convergence_iteration(starts(iC),0.99);
            ExtendedK99(q)=convergence_iteration(starts(iE),0.99);
            CandidateRestartStd(q)=std(vals(1:o.CandidateStarts),0,'omitnan');
            ExtendedRestartStd(q)=std(vals(1:o.ExtendedStarts),0,'omitnan');
            CandidateMaxFVGap(q)=max(abs([starts(1:o.CandidateStarts).validationGap]));
            ExtendedMaxFVGap(q)=max(abs([starts(1:o.ExtendedStarts).validationGap]));
            CandidateSelectionChanged(q)=iFastC~=iC; ExtendedSelectionChanged(q)=iFastE~=iE;
            CandidateSerialSeconds(q)=sum([starts(1:o.CandidateStarts).runtime])+ ...
                sum([starts(1:o.CandidateStarts).validationRuntime]);
            ExtendedSerialSeconds(q)=sum([starts(1:o.ExtendedStarts).runtime])+ ...
                sum([starts(1:o.ExtendedStarts).validationRuntime]);
            ExtendedOverCandidateRuntimeRatio(q)=ExtendedSerialSeconds(q)/max(CandidateSerialSeconds(q),eps);
            MeasuredExtendedWallSeconds(q)=c.runtime.wallSeconds;
            MaxAbsConstellationDelta(q)=max(abs(starts(iE).x_best(:)-starts(iC).x_best(:)));
            CaseFile(q)=string(path);
        end
    end

    T=table(SigmaX2,Replicate,Seed,SelectedDt,IsRefined,RuleVersion, ...
        BestCandidateAMI,BestExtendedAMI,ExtraValidatedGain,BestCandidateStart,BestExtendedStart, ...
        ExtraStartWins,MaterialExtraGain,CandidateK99,ExtendedK99,CandidateRestartStd,ExtendedRestartStd, ...
        CandidateMaxFVGap,ExtendedMaxFVGap,CandidateSelectionChanged,ExtendedSelectionChanged, ...
        CandidateSerialSeconds,ExtendedSerialSeconds,ExtendedOverCandidateRuntimeRatio, ...
        MeasuredExtendedWallSeconds,MaxAbsConstellationDelta,CaseFile);

    bySigma=aggregate_by_sigma(T,sigVec);
    summary=struct();
    summary.candidateStarts=double(o.CandidateStarts);
    summary.extendedStarts=double(o.ExtendedStarts);
    summary.nRows=height(T);
    summary.extraRestartTolerance=double(o.ExtraRestartTolerance);
    summary.fastValidationTolerance=double(o.FastValidationTolerance);
    summary.maxExtraValidatedGain=max(T.ExtraValidatedGain,[],'omitnan');
    summary.meanExtraValidatedGain=mean(T.ExtraValidatedGain,'omitnan');
    summary.nExtraStartWins=sum(T.ExtraStartWins);
    summary.nMaterialExtraGains=sum(T.MaterialExtraGain);
    summary.maxCandidateFVGap=max(T.CandidateMaxFVGap,[],'omitnan');
    summary.maxExtendedFVGap=max(T.ExtendedMaxFVGap,[],'omitnan');
    summary.nCandidateSelectionChanges=sum(T.CandidateSelectionChanged);
    summary.nExtendedSelectionChanges=sum(T.ExtendedSelectionChanged);
    summary.meanExtendedOverCandidateRuntimeRatio=mean(T.ExtendedOverCandidateRuntimeRatio,'omitnan');
    summary.maxCandidateK99=max(T.CandidateK99,[],'omitnan');
    summary.maxExtendedK99=max(T.ExtendedK99,[],'omitnan');
    summary.allRefined=all(T.IsRefined & abs(T.SelectedDt-0.005)<1e-14 & T.RuleVersion=="I10.3-v1");
    summary.extraRestartPass=summary.maxExtraValidatedGain<=o.ExtraRestartTolerance;
    summary.fastValidationPass=summary.maxExtendedFVGap<=o.FastValidationTolerance;
    summary.overallPass=summary.allRefined && summary.extraRestartPass && summary.fastValidationPass;
    summary.bySigma=bySigma;

    d=struct(); d.meta=fso_result_utils.run_metadata(mfilename);
    d.run=struct('id',runId,'outputDirectory',outDir,'productionRunDirectory',b.run.outputDirectory);
    d.settings=struct('M',double(o.M),'SNRdB',double(o.SNRdB),'TurbulenceVec',sigVec, ...
        'candidateStarts',double(o.CandidateStarts),'extendedStarts',double(o.ExtendedStarts), ...
        'nReplicates',nRep,'saMaxIter',double(o.saMaxIter),'T0',double(o.T0),'Tf',double(o.Tf), ...
        'BaseStd0',double(o.BaseStd0),'ItersPerTemp',double(o.ItersPerTemp), ...
        'ExtraRestartTolerance',double(o.ExtraRestartTolerance),'FastValidationTolerance',double(o.FastValidationTolerance));
    d.results=T; d.summary=summary; d.productionBundle=b;

    save(fullfile(outDir,'i10p4M32RestartReliability.mat'),'d','-v7.3');
    if o.WriteCSV
        writetable(T,fullfile(outDir,'i10p4M32RestartReliability_raw.csv'));
        writetable(bySigma,fullfile(outDir,'i10p4M32RestartReliability_bySigma.csv'));
    end

    fprintf('\nI.10.4 M32 RESTART RELIABILITY SUMMARY\n');
    fprintf('candidate/extended starts               : %d / %d\n',o.CandidateStarts,o.ExtendedStarts);
    fprintf('rows                                     : %d\n',height(T));
    fprintf('max/mean extra validated gain (6-4)      : %.3e / %.3e bit/symbol\n', ...
        summary.maxExtraValidatedGain,summary.meanExtraValidatedGain);
    fprintf('extra-start winners / material gains     : %d / %d\n',summary.nExtraStartWins,summary.nMaterialExtraGains);
    fprintf('max candidate/extended F-V gap           : %.3e / %.3e\n',summary.maxCandidateFVGap,summary.maxExtendedFVGap);
    fprintf('candidate/extended selection changes     : %d / %d\n',summary.nCandidateSelectionChanges,summary.nExtendedSelectionChanges);
    fprintf('mean extended/candidate serial-time ratio : %.3fx\n',summary.meanExtendedOverCandidateRuntimeRatio);
    fprintf('max candidate/extended k99               : %.0f / %.0f of %d\n', ...
        summary.maxCandidateK99,summary.maxExtendedK99,o.saMaxIter);
    fprintf('I.10.3 refined-policy metadata pass       : %d\n',summary.allRefined);
    fprintf('4-start extra-restart criterion pass      : %d\n',summary.extraRestartPass);
    fprintf('OVERALL I.10.4 PASS                       : %d\n',summary.overallPass);
    fprintf('\nRepeatability by turbulence (report-only):\n'); disp(bySigma);
end

function A=aggregate_by_sigma(T,sigVec)
    n=numel(sigVec); SigmaX2=nan(n,1); N=nan(n,1); CandidateMean=nan(n,1); CandidateStd=nan(n,1);
    CandidateRange=nan(n,1); ExtendedMean=nan(n,1); ExtendedStd=nan(n,1); ExtendedRange=nan(n,1);
    MaxExtraGain=nan(n,1); NExtraWins=nan(n,1); NMaterialExtra=nan(n,1); MeanRuntimeRatio=nan(n,1);
    for k=1:n
        R=T(abs(T.SigmaX2-sigVec(k))<1e-12,:); c=R.BestCandidateAMI; e=R.BestExtendedAMI;
        SigmaX2(k)=sigVec(k); N(k)=height(R); CandidateMean(k)=mean(c,'omitnan'); CandidateStd(k)=std(c,0,'omitnan');
        CandidateRange(k)=max(c)-min(c); ExtendedMean(k)=mean(e,'omitnan'); ExtendedStd(k)=std(e,0,'omitnan');
        ExtendedRange(k)=max(e)-min(e); MaxExtraGain(k)=max(R.ExtraValidatedGain,[],'omitnan');
        NExtraWins(k)=sum(R.ExtraStartWins); NMaterialExtra(k)=sum(R.MaterialExtraGain);
        MeanRuntimeRatio(k)=mean(R.ExtendedOverCandidateRuntimeRatio,'omitnan');
    end
    A=table(SigmaX2,N,CandidateMean,CandidateStd,CandidateRange,ExtendedMean,ExtendedStd,ExtendedRange, ...
        MaxExtraGain,NExtraWins,NMaterialExtra,MeanRuntimeRatio);
end

function k=convergence_iteration(s,fraction)
    improvement=s.bestMIFast-s.initialMI;
    if ~(isfinite(improvement)&&improvement>0), k=0; return; end
    target=s.initialMI+fraction*improvement; idx=find(s.history.bestMI>=target,1,'first');
    if isempty(idx), k=s.history.iter(end); else, k=s.history.iter(idx); end
end
function tf=positive_scalar(v),tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&isreal(v)&&v>0;end
function tf=nonnegative_scalar(v),tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&isreal(v)&&v>=0;end
function tf=finite_scalar(v),tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&isreal(v);end
function tf=positive_integer(v),tf=positive_scalar(v)&&mod(v,1)==0;end
function tf=nonnegative_vector(v),tf=isnumeric(v)&&isvector(v)&&~isempty(v)&&all(isfinite(v))&&all(v>=0);end
function t=sanitize_label(v),t=strtrim(char(string(v)));if isempty(t),t='run';return;end,t=regexprep(t,'[^A-Za-z0-9._-]+','_');t=regexprep(t,'_+','_');t=regexprep(t,'^[._-]+|[._-]+$','');if isempty(t),t='run';end,end
function [d,id]=make_unique_dir(parent,id),if ~exist(parent,'dir'),mkdir(parent);end,d=fullfile(parent,id);base=id;n=1;while exist(d,'dir'),n=n+1;id=sprintf('%s_%02d',base,n);d=fullfile(parent,id);end,[ok,msg]=mkdir(d);if ~ok,error('sim_i10p4_m32_restart_reliability:mkdir','%s',msg);end,end
function token=num_token(v,nDec),token=sprintf(['%0.' num2str(nDec) 'f'],double(v));token=strrep(token,'-','m');token=strrep(token,'+','p');token=strrep(token,'.','p');end
