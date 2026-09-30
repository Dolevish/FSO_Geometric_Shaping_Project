function cfg = i11_production_config()
%I11_PRODUCTION_CONFIG  Scientific freeze for the final IEEE-revision run.
%
% I.11 freezes the numerical/optimization design.  Operational controls
% (worker count, output directory, resume path) are intentionally excluded
% from the scientific signature and may change without changing the result.

    cfg=struct();
    cfg.freezeVersion='I11-v1';
    cfg.description='Final production freeze: flat restart tasks + checkpoint/resume';

    cfg.MVec=[4 8 16 32];
    cfg.P_avg=1;
    cfg.SNRVec=5:5:30;
    cfg.TurbulenceVec=[0 0.1 0.2 0.3];

    cfg.minGap=0.01;
    cfg.pinZero=true;

    cfg.FastFadingMethod='logh';
    cfg.LoghDtPolicy='case-adaptive';
    cfg.loghDt=0.01;
    cfg.loghDtRefined=0.005;
    cfg.loghSpanSigma=8;
    cfg.loghYBlockSize=512;
    cfg.YGridMode='candidate-adaptive';
    cfg.YGridScale=1.1;
    cfg.dtRuleVersion='I10.3-v1';
    cfg.runtimeTuningVersion='I10.2.2-v1';

    cfg.saMaxIter=4000;
    cfg.RestartPolicy='i11';
    cfg.RestartsByM=[6 6 8 8];
    cfg.nReplicates=2;
    cfg.T0=0.4;
    cfg.Tf=1e-3;
    cfg.BaseStd0=0.15;
    cfg.ItersPerTemp=50;
    cfg.baseSeed=20260928;
    cfg.historyEvery=50;

    cfg.parallelGranularity='restart-task';
    cfg.longestFirst=true;
    cfg.recommendedWorkers=8;
    cfg.workerNote=['8 workers was the fastest tested count on the local M=32 ' ...
        '200-iteration scaling diagnostic; worker count is operational, not scientific.'];

    nPhysical=numel(cfg.MVec)*numel(cfg.SNRVec)*numel(cfg.TurbulenceVec);
    cfg.nPhysicalCases=nPhysical;
    cfg.nOptimizationCases=nPhysical*cfg.nReplicates;
    cfg.totalRestartTasks=0;
    for k=1:numel(cfg.MVec)
        cfg.totalRestartTasks=cfg.totalRestartTasks + ...
            numel(cfg.SNRVec)*numel(cfg.TurbulenceVec)*cfg.nReplicates*cfg.RestartsByM(k);
    end

    cfg.restartDecisionNote=['Six restarts are used for M=4,8 and eight for M=16,32. ' ...
        'This is a conservative engineering choice. I.10.4 established that four starts ' ...
        'can be insufficient for tested M=32 cases; the I.11 policy is not claimed to be ' ...
        'the minimum sufficient restart count.'];
    cfg.globalOptimumClaim=false;
end
