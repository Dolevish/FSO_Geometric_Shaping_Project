function cfg = i11p1_production_config()
%I11P1_PRODUCTION_CONFIG  Revised frozen profile with an M-dependent SA budget.
%
% This supersedes I.11-v1 only for the canonical final production run. The
% earlier I.11 configuration remains in the repository as a historical
% record and continues to support its regression tests.

    cfg=i11_production_config();
    cfg.freezeVersion='I11.1-v1';
    cfg.description=['Final production freeze: checkpointed flat restarts, ' ...
        '6/6/8/8 restart map, and M-dependent iteration budget'];

    if isfield(cfg,'saMaxIter')
        cfg=rmfield(cfg,'saMaxIter');
    end

    cfg.IterationPolicy='i11p1';
    cfg.MaxIterByM=zeros(size(cfg.MVec));
    cfg.IterationPolicyMeta=cell(size(cfg.MVec));
    for k=1:numel(cfg.MVec)
        [cfg.MaxIterByM(k),cfg.IterationPolicyMeta{k}]=i11p1_iteration_policy(cfg.MVec(k));
    end
    cfg.iterationPolicyVersion='I11.1-v1';

    casesPerM=numel(cfg.SNRVec)*numel(cfg.TurbulenceVec)*cfg.nReplicates;
    cfg.totalSAIterations=0;
    for k=1:numel(cfg.MVec)
        cfg.totalSAIterations=cfg.totalSAIterations + ...
            casesPerM*cfg.RestartsByM(k)*cfg.MaxIterByM(k);
    end

    cfg.stageOrder={'m32','standard'};
    cfg.stageDecisionNote=['The canonical run executes the M=32 stage first so the ' ...
        'bottleneck always has the full flat restart-task worker pool. The remaining ' ...
        'M=4,8,16 cases are then executed as a second flat stage.'];
    cfg.iterationDecisionNote=['M=32 uses 8000 iterations/restart; M=4,8,16 use 4000. ' ...
        'The M=32 extension is a conservative production budget choice and is not ' ...
        'claimed to be the minimum sufficient or globally optimal budget.'];
    cfg.globalOptimumClaim=false;
end
