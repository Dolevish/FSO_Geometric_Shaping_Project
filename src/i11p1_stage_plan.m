function stages = i11p1_stage_plan(cfg)
%I11P1_STAGE_PLAN  Two-stage execution plan for the I.11.1 production freeze.
%
% Stage 1: M=32 only, 8000 iterations/restart.
% Stage 2: M=4,8,16, 4000 iterations/restart.
%
% Each stage still uses the validated single-level flat restart scheduler.
% Running M=32 first keeps the bottleneck saturated across the process pool
% instead of leaving a serial M=32 tail after lighter cases finish.

    if nargin<1 || isempty(cfg), cfg=i11p1_production_config(); end

    stages=repmat(struct('name','','MVec',[],'saMaxIter',NaN, ...
        'nOptimizationCases',NaN,'nRestartTasks',NaN,'totalSAIterations',NaN),1,2);

    stages(1)=make_stage('m32',32,8000,cfg);
    stages(2)=make_stage('standard',[4 8 16],4000,cfg);
end

function s=make_stage(name,MVec,maxIter,cfg)
    mask=ismember(cfg.MVec,MVec);
    if ~all(ismember(MVec,cfg.MVec))
        error('i11p1_stage_plan:UnknownM','Stage contains M outside the frozen production grid.');
    end
    casesPerM=numel(cfg.SNRVec)*numel(cfg.TurbulenceVec)*cfg.nReplicates;
    starts=cfg.RestartsByM(mask);

    s=struct();
    s.name=name;
    s.MVec=MVec;
    s.saMaxIter=maxIter;
    s.nOptimizationCases=casesPerM*numel(MVec);
    s.nRestartTasks=casesPerM*sum(starts);
    s.totalSAIterations=casesPerM*sum(starts)*maxIter;
end
