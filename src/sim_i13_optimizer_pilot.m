function result = sim_i13_optimizer_pilot(varargin)
%SIM_I13_OPTIMIZER_PILOT  Reviewer-2 SA-vs-Pattern-Search pilot.
%
% Runs Pattern Search on representative frozen I.11.1 channel cases and
% compares its validated winner against the corresponding production SA
% replicate. The raw restart initializations and scientific AMI evaluators
% match the SA setup.
%
% This is a PILOT: default BudgetFraction=0.05 is for runtime/behavior
% characterization only. Reviewer-facing conclusions require a later frozen
% benchmark with BudgetFraction=1 and a pre-declared case grid/tie tolerance.

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'ProductionBundlePath','',@text_scalar);
    addParameter(p,'ResultsRoot',fso_result_utils.default_results_root(),@text_scalar);
    addParameter(p,'OutputDirectory','',@text_scalar);
    addParameter(p,'BudgetFraction',0.05,@fraction_scalar);
    addParameter(p,'Replicate',1,@(v)isnumeric(v)&&isscalar(v)&&ismember(v,[1 2]));
    addParameter(p,'Cases',[8 20 0.1;16 20 0.2;32 20 0.3;32 30 0.2],@case_matrix);
    parse(p,varargin{:}); o=p.Results;

    A=analyze_i11p1_production_results( ...
        'ProductionBundlePath',o.ProductionBundlePath, ...
        'ResultsRoot',o.ResultsRoot,'WriteCSV',true);
    S=load(A.sourceBundlePath,'bundle');
    bProd=S.bundle;
    F=i11p1_production_config();

    if strlength(string(o.OutputDirectory))>0
        outDir=char(string(o.OutputDirectory));
    else
        outDir=fullfile(fileparts(A.outputDirectory),'i13_reviewer12','optimizer_pilot');
    end
    if exist(outDir,'dir')~=7,mkdir(outDir);end

    cases=double(o.Cases);
    n=size(cases,1);
    rows=cell(n,1);
    details=cell(n,1);

    fprintf('\n============================================================\n');
    fprintf('I.13C REVIEWER-2 OPTIMIZER PILOT\n');
    fprintf('Pattern Search vs frozen SA replicate %d\n',o.Replicate);
    fprintf('BudgetFraction=%.3f (pilot only; not reviewer-final)\n',o.BudgetFraction);
    fprintf('============================================================\n');

    for k=1:n
        M=cases(k,1); snr=cases(k,2); sig=cases(k,3);
        R=A.caseTable(A.caseTable.M==M & abs(A.caseTable.SNRdB-snr)<1e-12 & ...
            abs(A.caseTable.SigmaX2-sig)<1e-12 & A.caseTable.Replicate==o.Replicate,:);
        if height(R)~=1
            error('sim_i13_optimizer_pilot:CaseLookup', ...
                'Expected one production row for M=%g SNR=%g sigma=%g rep=%d.', ...
                M,snr,sig,o.Replicate);
        end

        [xBaseline,baselineAMI]=production_baseline(bProd,M,snr,sig);
        [maxIter,~]=i11p1_iteration_policy(M);
        nStarts=R.NStarts;
        fullPerStart=maxIter+1; % initial score + one proposal per SA iteration
        pilotPerStart=max(20,ceil(o.BudgetFraction*fullPerStart));

        cfg=build_fso_config(M,F.P_avg,snr,sig, ...
            'FastFadingMethod','logh', ...
            'loghDt',R.SelectedDt, ...
            'loghSpanSigma',F.loghSpanSigma, ...
            'loghYBlockSize',F.loghYBlockSize, ...
            'YGridMode','candidate-adaptive', ...
            'YGridScale',F.YGridScale, ...
            'minGap',F.minGap,'pinZero',F.pinZero, ...
            'saMaxIter',maxIter,'saNStarts',nStarts, ...
            'seedInit',R.Seed,'logEvery',0,'historyEvery',F.historyEvery);

        fprintf('\nCase %d/%d: M=%d SNR=%.1f sigma_R^2=%.1f | starts=%d | PS eval/start=%d (full fair=%d)\n', ...
            k,n,M,snr,sig,nStarts,pilotPerStart,fullPerStart);

        [psOut,psRuns]=patternsearch_multistart_validated( ...
            cfg,xBaseline,R.Seed,nStarts,pilotPerStart);

        row=struct();
        row.M=M; row.SNRdB=snr; row.SigmaR2=sig; row.Replicate=o.Replicate;
        row.Seed=R.Seed; row.NStarts=nStarts;
        row.SAIterationsPerStart=maxIter;
        row.FullFairEvalsPerStart=fullPerStart;
        row.PSPilotEvalsPerStart=pilotPerStart;
        row.BudgetFraction=o.BudgetFraction;
        row.PAMValidated=baselineAMI;
        row.SAValidated=R.AMIValidated;
        row.PSValidated=psOut.bestMIValidated;
        row.SAMinusPS=R.AMIValidated-psOut.bestMIValidated;
        row.SAGainOverPAM=R.AMIValidated-baselineAMI;
        row.PSGainOverPAM=psOut.bestMIValidated-baselineAMI;
        row.PSTotalFunctionEvaluations=psOut.totalFunctionEvaluations;
        row.PSTotalRuntimeSeconds=psOut.totalRuntime;
        row.PSMaxFastValidatorGap=psOut.maxAbsFastValidationGap;
        rows{k}=row;
        details{k}=struct('case',row,'patternsearch',psOut,'runs',psRuns);
    end

    summary=struct2table(vertcat(rows{:}));
    result=struct();
    result.version='I13C-pilot-v1';
    result.isReviewerFinal=false;
    result.reasonNotFinal='BudgetFraction is below 1 unless explicitly overridden; final case grid/tie tolerance not frozen yet.';
    result.sourceProduction=A.sourceBundlePath;
    result.summary=summary;
    result.details=details;
    result.outputDirectory=outDir;

    save(fullfile(outDir,'i13c_optimizer_pilot.mat'),'result','-v7.3');
    writetable(summary,fullfile(outDir,'i13c_optimizer_pilot.csv'));

    fprintf('\nI.13C PILOT SUMMARY (not reviewer-final)\n');
    disp(summary(:,{'M','SNRdB','SigmaR2','SAValidated','PSValidated','SAMinusPS', ...
        'PSTotalFunctionEvaluations','PSTotalRuntimeSeconds'}));
    fprintf('Outputs: %s\n',outDir);
end

function [x,ami]=production_baseline(b,M,snr,sig)
    if M==32,stage=b.stages.m32;else,stage=b.stages.standard;end
    iM=find(stage.axes.M==M,1);
    iS=find(abs(stage.axes.SNR_dB-snr)<1e-12,1);
    iT=find(abs(stage.axes.sigma_X_sq-sig)<1e-12,1);
    if isempty(iM)||isempty(iS)||isempty(iT)
        error('sim_i13_optimizer_pilot:BaselineLookup','Production baseline not found.');
    end
    z=stage.baseline{iM,iT,iS};
    x=z.x(:);
    if isfield(z,'amiValidated')
        ami=z.amiValidated;
    elseif isfield(z,'AMIValidated')
        ami=z.AMIValidated;
    else
        error('sim_i13_optimizer_pilot:BaselineAMI','Validated baseline AMI field not found.');
    end
end

function tf=text_scalar(v),tf=ischar(v)||(isstring(v)&&isscalar(v));end
function tf=fraction_scalar(v),tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>0&&v<=1;end
function tf=case_matrix(v),tf=isnumeric(v)&&ismatrix(v)&&size(v,2)==3&&all(isfinite(v(:)));end
