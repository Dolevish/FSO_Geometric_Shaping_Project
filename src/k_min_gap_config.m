function K = k_min_gap_config()
%K_MIN_GAP_CONFIG  Frozen Commit-K minimum-spacing sensitivity protocol.
%
% Scientific grid:
%   M            = [16 32]
%   SNR          = 20 dB
%   sigma_R^2    = 0.2
%   d_min        = [0 0.005 0.01 0.02 0.05]
%   replicates   = 2
%   restarts     = 8 per replicate
%   iterations   = 8000 per restart for BOTH M=16 and M=32
%
% Commit K deliberately gives M=16 the same 8000-iteration budget as M=32
% so d_min sensitivity is compared under a symmetric conservative budget.
% All other numerical/SA settings are inherited from the frozen I.11.1
% production profile. No global-optimum claim is implied.

    P=i11p1_production_config();

    K=struct();
    K.version='K-v1';
    K.description='Canonical minimum-spacing sensitivity rerun using frozen production numerics';
    K.MVec=[16 32];
    K.SNRdB=20;
    K.SigmaR2=0.2;
    K.MinGapVec=[0 0.005 0.01 0.02 0.05];
    K.Replicates=[1 2];
    K.NStarts=8;
    K.MaxIter=8000;
    K.P_avg=P.P_avg;
    K.pinZero=P.pinZero;
    K.baseSeed=P.baseSeed;

    % Frozen production evaluator.
    K.FastFadingMethod=P.FastFadingMethod;
    K.LoghDtPolicy=P.LoghDtPolicy;
    K.loghDt=P.loghDt;
    K.loghDtRefined=P.loghDtRefined;
    K.loghSpanSigma=P.loghSpanSigma;
    K.loghYBlockSize=P.loghYBlockSize;
    K.YGridMode=P.YGridMode;
    K.YGridScale=P.YGridScale;
    K.historyEvery=P.historyEvery;

    % Frozen production SA profile.
    K.T0=P.T0;
    K.Tf=P.Tf;
    K.BaseStd0=P.BaseStd0;
    K.ItersPerTemp=P.ItersPerTemp;

    % Commit-K parallel execution policy.
    K.ParallelWorkers=8;
    K.scheduler='parfeval-dynamic-v1';
    K.taskGranularity='one validated SA restart';
    K.longestFirst=true;
    K.noNestedParallelism=true;

    % Verify and freeze the I.10.3 log-h resolution for the two cases.
    K.SelectedDtByM=zeros(size(K.MVec));
    K.DtMetaByM=cell(size(K.MVec));
    for i=1:numel(K.MVec)
        [K.SelectedDtByM(i),K.DtMetaByM{i}]=select_logh_dt( ...
            K.MVec(i),K.SNRdB,K.SigmaR2, ...
            'Policy',K.LoghDtPolicy,'BaseDt',K.loghDt, ...
            'RefinedDt',K.loghDtRefined,'FixedDt',K.loghDt);
    end
    if any(abs(K.SelectedDtByM-0.01)>1e-12)
        error('k_min_gap_config:UnexpectedDt', ...
            'Commit K expects dt=0.01 for both SNR=20 dB cases.');
    end

    K.nPhysicalSensitivityPoints=numel(K.MVec)*numel(K.MinGapVec);
    K.nReplicateCases=K.nPhysicalSensitivityPoints*numel(K.Replicates);
    K.nRestartTasks=K.nReplicateCases*K.NStarts;
    K.totalSAIterations=K.nRestartTasks*K.MaxIter;
    K.nominalFastObjectiveEvaluations=K.nRestartTasks*(K.MaxIter+1);
    K.globalOptimumClaim=false;
end
