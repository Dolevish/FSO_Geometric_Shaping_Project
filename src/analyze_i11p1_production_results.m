function analysis = analyze_i11p1_production_results(varargin)
%ANALYZE_I11P1_PRODUCTION_RESULTS  Post-production analysis for I.11.1.
%
%   analysis = analyze_i11p1_production_results()
%   analysis = analyze_i11p1_production_results('ProductionBundlePath',path)
%
% This is the canonical I.12A post-processing entry point for the frozen
% I.11.1 production run. It does not recompute AMI and does not run SA.
% Every number is read from i11p1ProductionBundle.mat.
%
% Outputs:
%   i12_physical_summary.csv
%       One row per (M,sigma_X^2,SNR), aggregating the two replicates.
%   i12_per_M_summary.csv
%       Compact repeatability / accuracy / gain summary by M.
%   i12_tableI_M8_SNR20.csv
%       Candidate replacement for manuscript Table I.
%   i12_selection_changes.csv
%       Optimization cases where fast-winner and validated-winner differ.
%   i12_repeatability_sorted.csv
%       Physical points sorted by replicate-to-replicate AMI range.
%   i12_fig4_M8_SNR20_constellations.csv
%       PAM, both GS replicates, and a deterministic representative
%       constellation for the four turbulence values at M=8,SNR=20 dB.
%   i12_analysis.mat
%       Complete analysis structure.
%
% Replicate aggregation:
%   Performance curves/tables are summarized by the arithmetic mean of the
%   two independently validated replicate winners, with sample standard
%   deviation retained separately. No "best replicate" aggregation is used.
%
% Representative constellation policy:
%   For geometry-only inspection, the representative replicate is the one
%   whose validated AMI is closest to the two-replicate mean; ties select
%   the lower replicate index. This avoids cherry-picking the best replicate
%   and is diagnostic only until the paper-presentation policy is frozen.
%
% No global-optimum claim is implied.

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'ProductionBundlePath','',@text_scalar);
    addParameter(p,'ResultsRoot',fso_result_utils.default_results_root(),@text_scalar);
    addParameter(p,'OutputDirectory','',@text_scalar);
    addParameter(p,'WriteCSV',true,@logical_scalar);
    parse(p,varargin{:}); o=p.Results;

    bundlePath=resolve_bundle_path(o.ProductionBundlePath,o.ResultsRoot);
    S=load(bundlePath,'bundle');
    if ~isfield(S,'bundle')
        error('analyze_i11p1_production_results:MissingBundle', ...
            'File does not contain variable bundle: %s',bundlePath);
    end
    b=S.bundle;
    F=i11p1_production_config();
    validate_bundle(b,F,bundlePath);

    if strlength(string(o.OutputDirectory))>0
        outDir=char(string(o.OutputDirectory));
    else
        outDir=fullfile(fileparts(bundlePath),'i12_analysis');
    end
    if exist(outDir,'dir')~=7
        [ok,msg]=mkdir(outDir);
        if ~ok,error('analyze_i11p1_production_results:mkdir','%s',msg);end
    end

    T=b.summary.caseTable;
    T=sortrows(T,{'M','SigmaX2','SNRdB','Replicate'});

    physical=build_physical_summary(T,F);
    perM=build_per_m_summary(physical,T,F);
    tableI=physical(physical.M==8 & abs(physical.SNRdB-20)<1e-12, ...
        {'SigmaX2','PAMValidated','GSMeanValidated','GSStdAcrossReplicates', ...
         'GainMeanBits','GainStdAcrossReplicates','ReplicateRangeBits'});
    tableI=sortrows(tableI,'SigmaX2');

    selectionChanges=T(logical(T.SelectionChanged),:);
    repeatabilitySorted=sortrows(physical, ...
        {'ReplicateRangeBits','M','SigmaX2','SNRdB'}, ...
        {'descend','ascend','ascend','ascend'});

    fig4Table=build_m8_snr20_constellation_table(b,physical);

    [maxGap,maxGapRow]=max(T.MaxAbsFastValidationGap);
    maxGapCase=T(maxGapRow,{'M','SigmaX2','SNRdB','Replicate', ...
        'MaxAbsFastValidationGap','SelectionChanged'});
    [maxRepRange,maxRepRow]=max(physical.ReplicateRangeBits);
    maxRepCase=physical(maxRepRow,{'M','SigmaX2','SNRdB', ...
        'GSMeanValidated','GSStdAcrossReplicates','ReplicateRangeBits'});

    analysis=struct();
    analysis.version='I12A-v1';
    analysis.sourceBundlePath=bundlePath;
    analysis.outputDirectory=outDir;
    analysis.productionFreeze=F;
    analysis.aggregationPolicy='mean validated AMI across replicates; sample std retained';
    analysis.representativeConstellationPolicy=[ ...
        'replicate whose validated AMI is closest to the across-replicate mean; ' ...
        'tie -> lower replicate index; diagnostic only'];
    analysis.caseTable=T;
    analysis.physicalTable=physical;
    analysis.perMTable=perM;
    analysis.tableI_M8_SNR20=tableI;
    analysis.selectionChanges=selectionChanges;
    analysis.repeatabilitySorted=repeatabilitySorted;
    analysis.fig4_M8_SNR20=fig4Table;
    analysis.metrics=struct( ...
        'nOptimizationCases',height(T), ...
        'nPhysicalPoints',height(physical), ...
        'nSelectionChanges',height(selectionChanges), ...
        'selectionChangeFraction',height(selectionChanges)/height(T), ...
        'maxFastValidatorGap',maxGap, ...
        'maxFastValidatorGapCase',maxGapCase, ...
        'maxReplicateRangeBits',maxRepRange, ...
        'maxReplicateRangeCase',maxRepCase, ...
        'minIndividualGainBits',min(T.GainBits), ...
        'minMeanPhysicalGainBits',min(physical.GainMeanBits), ...
        'nNegativeIndividualGains',sum(T.GainBits<0), ...
        'nNegativeMeanPhysicalGains',sum(physical.GainMeanBits<0));

    save(fullfile(outDir,'i12_analysis.mat'),'analysis','-v7.3');
    if o.WriteCSV
        writetable(physical,fullfile(outDir,'i12_physical_summary.csv'));
        writetable(perM,fullfile(outDir,'i12_per_M_summary.csv'));
        writetable(tableI,fullfile(outDir,'i12_tableI_M8_SNR20.csv'));
        writetable(selectionChanges,fullfile(outDir,'i12_selection_changes.csv'));
        writetable(repeatabilitySorted,fullfile(outDir,'i12_repeatability_sorted.csv'));
        writetable(fig4Table,fullfile(outDir,'i12_fig4_M8_SNR20_constellations.csv'));
    end

    print_report(analysis);
end

function physical=build_physical_summary(T,F)
    keys=unique(T(:,{'M','SigmaX2','SNRdB'}),'rows','sorted');
    keys=sortrows(keys,{'M','SigmaX2','SNRdB'});
    n=height(keys);

    PAMValidated=nan(n,1);
    GSRep1Validated=nan(n,1);
    GSRep2Validated=nan(n,1);
    GSMeanValidated=nan(n,1);
    GSStdAcrossReplicates=nan(n,1);
    GainRep1Bits=nan(n,1);
    GainRep2Bits=nan(n,1);
    GainMeanBits=nan(n,1);
    GainStdAcrossReplicates=nan(n,1);
    ReplicateRangeBits=nan(n,1);
    RepresentativeReplicate=nan(n,1);
    MaxAbsFastValidationGap=nan(n,1);
    SelectionChanges=zeros(n,1);
    MeanRestartValidatedStd=nan(n,1);
    MaxRestartValidatedStd=nan(n,1);
    MeanK99=nan(n,1);
    MaxK99=nan(n,1);
    SelectedDt=nan(n,1);
    NStarts=nan(n,1);

    for k=1:n
        idx=T.M==keys.M(k) & ...
            abs(T.SigmaX2-keys.SigmaX2(k))<1e-12 & ...
            abs(T.SNRdB-keys.SNRdB(k))<1e-12;
        R=sortrows(T(idx,:),'Replicate');
        if height(R)~=F.nReplicates || ~isequal(R.Replicate(:).',(1:F.nReplicates))
            error('analyze_i11p1_production_results:ReplicateCount', ...
                'Expected replicates 1:%d for M=%d sigma=%.3f SNR=%.1f.', ...
                F.nReplicates,keys.M(k),keys.SigmaX2(k),keys.SNRdB(k));
        end

        pamByRep=R.AMIValidated-R.GainBits;
        if max(pamByRep)-min(pamByRep)>1e-9
            error('analyze_i11p1_production_results:BaselineMismatch', ...
                'Replicate baselines disagree for M=%d sigma=%.3f SNR=%.1f.', ...
                keys.M(k),keys.SigmaX2(k),keys.SNRdB(k));
        end
        if max(R.SelectedDt)-min(R.SelectedDt)>1e-12
            error('analyze_i11p1_production_results:DtMismatch', ...
                'Replicate dt differs for one physical point.');
        end
        if numel(unique(R.NStarts))~=1
            error('analyze_i11p1_production_results:RestartMismatch', ...
                'Replicate restart counts differ for one physical point.');
        end

        gs=R.AMIValidated(:); gain=R.GainBits(:);
        gsMean=mean(gs);
        dRep=abs(gs-gsMean);
        repPos=find(dRep<=min(dRep)+1e-12,1,'first');

        PAMValidated(k)=mean(pamByRep);
        GSRep1Validated(k)=gs(1);
        GSRep2Validated(k)=gs(2);
        GSMeanValidated(k)=gsMean;
        GSStdAcrossReplicates(k)=std(gs,0);
        GainRep1Bits(k)=gain(1);
        GainRep2Bits(k)=gain(2);
        GainMeanBits(k)=mean(gain);
        GainStdAcrossReplicates(k)=std(gain,0);
        ReplicateRangeBits(k)=max(gs)-min(gs);
        RepresentativeReplicate(k)=R.Replicate(repPos);
        MaxAbsFastValidationGap(k)=max(R.MaxAbsFastValidationGap);
        SelectionChanges(k)=sum(logical(R.SelectionChanged));
        MeanRestartValidatedStd(k)=mean(R.RestartValidatedStd,'omitnan');
        MaxRestartValidatedStd(k)=max(R.RestartValidatedStd,[],'omitnan');
        MeanK99(k)=mean(R.K99,'omitnan');
        MaxK99(k)=max(R.K99,[],'omitnan');
        SelectedDt(k)=R.SelectedDt(1);
        NStarts(k)=R.NStarts(1);
    end

    physical=[keys table(PAMValidated,GSRep1Validated,GSRep2Validated, ...
        GSMeanValidated,GSStdAcrossReplicates,GainRep1Bits,GainRep2Bits, ...
        GainMeanBits,GainStdAcrossReplicates,ReplicateRangeBits, ...
        RepresentativeReplicate,MaxAbsFastValidationGap,SelectionChanges, ...
        MeanRestartValidatedStd,MaxRestartValidatedStd,MeanK99,MaxK99, ...
        SelectedDt,NStarts)];
end

function perM=build_per_m_summary(P,T,F)
    MVec=F.MVec(:); n=numel(MVec);
    MeanGainBits=nan(n,1); MinGainBits=nan(n,1); MaxGainBits=nan(n,1);
    MeanReplicateRangeBits=nan(n,1); MaxReplicateRangeBits=nan(n,1);
    MaxAbsFastValidationGap=nan(n,1); SelectionChanges=zeros(n,1);
    MeanRestartValidatedStd=nan(n,1); MaxRestartValidatedStd=nan(n,1);

    for k=1:n
        ip=P.M==MVec(k); it=T.M==MVec(k);
        MeanGainBits(k)=mean(P.GainMeanBits(ip));
        MinGainBits(k)=min(P.GainMeanBits(ip));
        MaxGainBits(k)=max(P.GainMeanBits(ip));
        MeanReplicateRangeBits(k)=mean(P.ReplicateRangeBits(ip));
        MaxReplicateRangeBits(k)=max(P.ReplicateRangeBits(ip));
        MaxAbsFastValidationGap(k)=max(T.MaxAbsFastValidationGap(it));
        SelectionChanges(k)=sum(logical(T.SelectionChanged(it)));
        MeanRestartValidatedStd(k)=mean(P.MeanRestartValidatedStd(ip),'omitnan');
        MaxRestartValidatedStd(k)=max(P.MaxRestartValidatedStd(ip),[],'omitnan');
    end
    M=MVec;
    perM=table(M,MeanGainBits,MinGainBits,MaxGainBits, ...
        MeanReplicateRangeBits,MaxReplicateRangeBits, ...
        MaxAbsFastValidationGap,SelectionChanges, ...
        MeanRestartValidatedStd,MaxRestartValidatedStd);
end

function C=build_m8_snr20_constellation_table(b,physical)
    stage=get_stage_for_m(b,8);
    iM=find(stage.axes.M==8,1);
    iSNR=find(abs(stage.axes.SNR_dB-20)<1e-12,1);
    sig=stage.axes.sigma_X_sq(:);
    nSig=numel(sig); M=8; nRep=numel(stage.axes.replicate);

    if isempty(iM)||isempty(iSNR)||nRep~=2
        error('analyze_i11p1_production_results:Fig4Axes', ...
            'Expected M=8, SNR=20 dB, and exactly two replicates.');
    end

    repChoice=nan(nSig,1);
    pam=nan(nSig,M); gs1=nan(nSig,M); gs2=nan(nSig,M); gsRep=nan(nSig,M);
    GSRep1AMI=nan(nSig,1); GSRep2AMI=nan(nSig,1); GSMeanAMI=nan(nSig,1);

    for iSig=1:nSig
        pidx=physical.M==8 & abs(physical.SigmaX2-sig(iSig))<1e-12 & ...
            abs(physical.SNRdB-20)<1e-12;
        if sum(pidx)~=1,error('analyze_i11p1_production_results:Fig4Physical','Missing M8/SNR20 physical row.');end
        repChoice(iSig)=physical.RepresentativeReplicate(pidx);
        GSRep1AMI(iSig)=physical.GSRep1Validated(pidx);
        GSRep2AMI(iSig)=physical.GSRep2Validated(pidx);
        GSMeanAMI(iSig)=physical.GSMeanValidated(pidx);

        base=stage.baseline{iM,iSig,iSNR};
        pam(iSig,:)=base.x(:).';
        gs1(iSig,:)=stage.summary.xWinner{iM,iSig,iSNR,1}(:).';
        gs2(iSig,:)=stage.summary.xWinner{iM,iSig,iSNR,2}(:).';
        gsRep(iSig,:)=stage.summary.xWinner{iM,iSig,iSNR,repChoice(iSig)}(:).';
    end

    C=table(sig,repChoice,GSRep1AMI,GSRep2AMI,GSMeanAMI, ...
        'VariableNames',{'SigmaX2','RepresentativeReplicate','GSRep1AMI','GSRep2AMI','GSMeanAMI'});
    for j=1:M
        C.(sprintf('PAM_x%d',j))=pam(:,j);
    end
    for j=1:M
        C.(sprintf('GS_rep1_x%d',j))=gs1(:,j);
    end
    for j=1:M
        C.(sprintf('GS_rep2_x%d',j))=gs2(:,j);
    end
    for j=1:M
        C.(sprintf('GS_representative_x%d',j))=gsRep(:,j);
    end
end

function stage=get_stage_for_m(b,M)
    if M==32
        stage=b.stages.m32;
    else
        stage=b.stages.standard;
    end
    if ~any(stage.axes.M==M)
        error('analyze_i11p1_production_results:StageAxes','M=%d is missing from stage.',M);
    end
end

function validate_bundle(b,F,path)
    need={'productionFreeze','stages','summary'};
    for k=1:numel(need)
        if ~isfield(b,need{k})
            error('analyze_i11p1_production_results:BadBundle', ...
                'Missing bundle.%s in %s.',need{k},path);
        end
    end
    if ~isfield(b.productionFreeze,'freezeVersion') || ...
            ~strcmp(b.productionFreeze.freezeVersion,F.freezeVersion)
        error('analyze_i11p1_production_results:FreezeVersion', ...
            'Bundle is not the frozen %s production profile.',F.freezeVersion);
    end
    if ~isfield(b.stages,'standard') || ~isfield(b.stages,'m32')
        error('analyze_i11p1_production_results:MissingStage','Both standard and m32 stages are required.');
    end
    if ~isfield(b.summary,'caseTable') || ~istable(b.summary.caseTable)
        error('analyze_i11p1_production_results:MissingCaseTable','Combined case table is missing.');
    end

    T=b.summary.caseTable;
    required={'M','SigmaX2','SNRdB','Replicate','NStarts','SelectedDt', ...
        'AMIValidated','GainBits','K99','RestartValidatedStd', ...
        'MaxAbsFastValidationGap','SelectionChanged'};
    missing=setdiff(required,T.Properties.VariableNames);
    if ~isempty(missing)
        error('analyze_i11p1_production_results:MissingColumns', ...
            'Case table is missing: %s',strjoin(missing,', '));
    end
    if height(T)~=F.nOptimizationCases
        error('analyze_i11p1_production_results:CaseCount', ...
            'Expected %d optimization cases, found %d.',F.nOptimizationCases,height(T));
    end
    if ~isequal(unique(T.M(:)).',F.MVec)
        error('analyze_i11p1_production_results:MGrid','Unexpected M grid.');
    end
    if ~isequal(unique(T.SNRdB(:)).',F.SNRVec)
        error('analyze_i11p1_production_results:SNRGrid','Unexpected SNR grid.');
    end
    if max(abs(unique(T.SigmaX2(:)).'-F.TurbulenceVec))>1e-12
        error('analyze_i11p1_production_results:TurbulenceGrid','Unexpected turbulence grid.');
    end
end

function path=resolve_bundle_path(requested,resultsRoot)
    if strlength(string(requested))>0
        path=char(string(requested));
        if exist(path,'file')~=2
            error('analyze_i11p1_production_results:BundleNotFound','Bundle not found: %s',path);
        end
        return;
    end

    base=char(string(resultsRoot));
    roots={fullfile(base,'i11p1_production'), ...
        fullfile(base,'ieee_revision','i11p1_production')};
    roots=unique(roots,'stable');

    D=[];
    searched=strings(0,1);
    for i=1:numel(roots)
        searched(end+1,1)=string(roots{i}); %#ok<AGROW>
        Di=dir(fullfile(roots{i},'*','i11p1ProductionBundle.mat'));
        if ~isempty(Di)
            if isempty(D)
                D=Di;
            else
                D=[D; Di]; %#ok<AGROW>
            end
        end
    end
    if isempty(D)
        error('analyze_i11p1_production_results:NoProductionBundle', ...
            ['No i11p1ProductionBundle.mat found. Searched:\n  %s\n' ...
             'Supply ''ProductionBundlePath'' explicitly if needed.'], ...
             strjoin(cellstr(searched),sprintf('\n  ')));
    end
    [~,k]=max([D.datenum]); %#ok<DATNM>
    path=fullfile(D(k).folder,D(k).name);
end

function print_report(A)
    P=A.physicalTable; M=A.metrics;
    fprintf('\n============================================================\n');
    fprintf('I.12A POST-PRODUCTION RESULTS ANALYSIS (%s)\n',A.version);
    fprintf('Source: %s\n',A.sourceBundlePath);
    fprintf('Cases: %d | physical points: %d | replicates: 2\n', ...
        M.nOptimizationCases,M.nPhysicalPoints);
    fprintf('Aggregation: mean validated AMI across replicates; sample std retained.\n');
    fprintf('Selection changes: %d/%d (%.2f%%)\n', ...
        M.nSelectionChanges,M.nOptimizationCases,100*M.selectionChangeFraction);
    fprintf('Max fast-validator gap: %.3e bits/symbol\n',M.maxFastValidatorGap);
    fprintf('Max replicate AMI range: %.6f bits/symbol\n',M.maxReplicateRangeBits);
    fprintf('Min individual gain: %.6f | min mean physical-point gain: %.6f bits/symbol\n', ...
        M.minIndividualGainBits,M.minMeanPhysicalGainBits);
    fprintf('Negative gains: %d individual cases | %d mean physical points\n', ...
        M.nNegativeIndividualGains,M.nNegativeMeanPhysicalGains);

    fprintf('\nPer-M summary:\n');
    disp(A.perMTable);

    fprintf('\nCandidate Table I: M=8, SNR=20 dB\n');
    disp(A.tableI_M8_SNR20);

    fprintf('\nWorst replicate-consistency point:\n');
    disp(M.maxReplicateRangeCase);
    fprintf('Worst fast-validator point:\n');
    disp(M.maxFastValidatorGapCase);

    fprintf('Outputs: %s\n',A.outputDirectory);
    fprintf('============================================================\n');
end

function tf=text_scalar(v),tf=ischar(v)||(isstring(v)&&isscalar(v));end
function tf=logical_scalar(v),tf=islogical(v)&&isscalar(v);end
