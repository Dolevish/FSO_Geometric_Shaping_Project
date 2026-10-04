function result = sim_i13_zhou_baseline_full_grid(varargin)
%SIM_I13_ZHOU_BASELINE_FULL_GRID  Reviewer-1 literature-GS comparison.
%
% Evaluates the deterministic Zhou-2025 optical-intensity GS construction
% with the SAME independent No-CSI FSO validator used for the frozen I.11.1
% production results.
%
% This function performs no SA and no channel-specific optimization.
% Therefore:
%   Zhou-PAM gain       = generic optical-intensity shaping contribution.
%   Proposed-Zhou gain  = additional gain from adapting geometry to the
%                         marginalized No-CSI FSO likelihood.
%
% Output: one row for each of the 96 frozen physical points, CSV + MAT.

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'ProductionBundlePath','',@text_scalar);
    addParameter(p,'ResultsRoot',fso_result_utils.default_results_root(),@text_scalar);
    addParameter(p,'OutputDirectory','',@text_scalar);
    addParameter(p,'ExtraBits',2,@(v)isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>=0&&mod(v,1)==0);
    parse(p,varargin{:}); o=p.Results;

    A=analyze_i11p1_production_results( ...
        'ProductionBundlePath',o.ProductionBundlePath, ...
        'ResultsRoot',o.ResultsRoot,'WriteCSV',true);

    if strlength(string(o.OutputDirectory))>0
        outDir=char(string(o.OutputDirectory));
    else
        outDir=fullfile(fileparts(A.outputDirectory),'i13_reviewer12','zhou_baseline');
    end
    if exist(outDir,'dir')~=7
        [ok,msg]=mkdir(outDir);
        if ~ok,error('sim_i13_zhou_baseline_full_grid:mkdir','%s',msg);end
    end

    F=i11p1_production_config();
    P=sortrows(A.physicalTable,{'M','SigmaX2','SNRdB'});
    if height(P)~=F.nPhysicalCases
        error('sim_i13_zhou_baseline_full_grid:PhysicalCount', ...
            'Expected %d physical points, found %d.',F.nPhysicalCases,height(P));
    end

    % Deterministic literature constellations: one per M, independent of
    % SNR/turbulence by construction.
    zhouByM=cell(numel(F.MVec),1);
    for iM=1:numel(F.MVec)
        zhouByM{iM}=zhou2025_optical_gs_constellation( ...
            F.MVec(iM),F.P_avg,'ExtraBits',o.ExtraBits,'Regularize',true);
    end

    n=height(P);
    ZhouValidated=nan(n,1);
    ZhouMinusPAM=nan(n,1);
    ProposedMinusZhou=nan(n,1);
    ProposedMinusPAM=P.GainMeanBits;
    ZhouValidationSeconds=nan(n,1);
    ZhouMinGap=nan(n,1);
    ZhouPAPR=nan(n,1);
    ZhouBits=nan(n,1);

    fprintf('\n============================================================\n');
    fprintf('I.13B ZHOU-2025 LITERATURE GS BASELINE\n');
    fprintf('Physical points: %d | validator-only, no optimization\n',n);
    fprintf('Regularization: b+n with n=%d extra bits\n',o.ExtraBits);
    fprintf('============================================================\n');

    for k=1:n
        M=P.M(k); snr=P.SNRdB(k); sig=P.SigmaX2(k);
        iM=find(F.MVec==M,1);
        Z=zhouByM{iM};

        % Recreate the frozen physical channel. Fast-evaluator parameters do
        % not affect AMI_Validator, but are supplied consistently.
        cfg=build_fso_config(M,F.P_avg,snr,sig, ...
            'FastFadingMethod','logh', ...
            'loghDt',P.SelectedDt(k), ...
            'loghSpanSigma',F.loghSpanSigma, ...
            'loghYBlockSize',F.loghYBlockSize, ...
            'YGridMode','candidate-adaptive', ...
            'YGridScale',F.YGridScale, ...
            'minGap',F.minGap, ...
            'pinZero',F.pinZero, ...
            'saMaxIter',1,'saNStarts',1,'logEvery',0);

        % Zhou n=2 is comfortably feasible for the frozen d_min=0.01 for
        % M=4,8,16,32. Fail loudly if this ever changes.
        AMI_functions.assert_constellation_feasible(Z.x,cfg,1e-10);

        t=tic;
        iz=cfg.AMI_Validator(Z.x);
        ZhouValidationSeconds(k)=toc(t);
        ZhouValidated(k)=iz;
        ZhouMinusPAM(k)=iz-P.PAMValidated(k);
        ProposedMinusZhou(k)=P.GSMeanValidated(k)-iz;
        ZhouMinGap(k)=Z.minGap;
        ZhouPAPR(k)=Z.PAPR;
        ZhouBits(k)=Z.bits;

        fprintf('%3d/%3d | M=%2d SNR=%2.0f sigma_R^2=%.1f | PAM=%.6f Zhou=%.6f Proposed=%.6f | Ggen=%+.6f Gadapt=%+.6f | %.2fs\n', ...
            k,n,M,snr,sig,P.PAMValidated(k),iz,P.GSMeanValidated(k), ...
            ZhouMinusPAM(k),ProposedMinusZhou(k),ZhouValidationSeconds(k));
    end

    comparison=[P(:,{'M','SigmaX2','SNRdB','PAMValidated','GSMeanValidated', ...
        'GSStdAcrossReplicates','GainMeanBits'}) ...
        table(ZhouValidated,ZhouMinusPAM,ProposedMinusZhou,ProposedMinusPAM, ...
              ZhouValidationSeconds,ZhouMinGap,ZhouPAPR,ZhouBits)];

    perM=groupsummary(comparison,'M',{'mean','min','max'}, ...
        {'ZhouMinusPAM','ProposedMinusZhou','ProposedMinusPAM'});

    % Focused Table-I operating point.
    tableI=comparison(comparison.M==8 & abs(comparison.SNRdB-20)<1e-12,:);
    tableI=sortrows(tableI,'SigmaX2');

    metrics=struct();
    metrics.nPhysicalPoints=height(comparison);
    metrics.meanGenericGain=mean(comparison.ZhouMinusPAM);
    metrics.meanAdaptationGain=mean(comparison.ProposedMinusZhou);
    metrics.minAdaptationGain=min(comparison.ProposedMinusZhou);
    metrics.maxAdaptationGain=max(comparison.ProposedMinusZhou);
    metrics.nProposedBelowZhou=sum(comparison.ProposedMinusZhou < -1e-10);
    metrics.nZhouBelowPAM=sum(comparison.ZhouMinusPAM < -1e-10);
    metrics.maxValidationSeconds=max(comparison.ZhouValidationSeconds);
    metrics.totalValidationSeconds=sum(comparison.ZhouValidationSeconds);

    result=struct();
    result.version='I13B-v1';
    result.sourceProduction=A.sourceBundlePath;
    result.extraBits=o.ExtraBits;
    result.zhouByM=zhouByM;
    result.comparison=comparison;
    result.perM=perM;
    result.tableI_M8_SNR20=tableI;
    result.metrics=metrics;
    result.outputDirectory=outDir;

    save(fullfile(outDir,'i13b_zhou_baseline.mat'),'result','-v7.3');
    writetable(comparison,fullfile(outDir,'i13b_zhou_full_grid.csv'));
    writetable(perM,fullfile(outDir,'i13b_zhou_per_M.csv'));
    writetable(tableI,fullfile(outDir,'i13b_zhou_M8_SNR20.csv'));

    fprintf('\nI.13B SUMMARY\n');
    fprintf('mean Zhou-PAM generic gain      : %+.6f bit/symbol\n',metrics.meanGenericGain);
    fprintf('mean Proposed-Zhou adapt. gain  : %+.6f bit/symbol\n',metrics.meanAdaptationGain);
    fprintf('min/max Proposed-Zhou gain      : %+.6f / %+.6f\n', ...
        metrics.minAdaptationGain,metrics.maxAdaptationGain);
    fprintf('Proposed below Zhou (>1e-10)    : %d/%d\n', ...
        metrics.nProposedBelowZhou,height(comparison));
    fprintf('Zhou below PAM (>1e-10)         : %d/%d\n', ...
        metrics.nZhouBelowPAM,height(comparison));
    fprintf('Outputs: %s\n',outDir);
end

function tf=text_scalar(v)
    tf=ischar(v)||(isstring(v)&&isscalar(v));
end
