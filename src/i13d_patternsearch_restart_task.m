function result = i13d_patternsearch_restart_task(cfg,xBaseline,masterSeed,restartIndex,maxEvals,psPolicy,identity,checkpointPath,signatureText)
%I13D_PATTERNSEARCH_RESTART_TASK  One checkpointed Pattern Search restart.
%
% Atomic task for the I.13D dynamic parfeval scheduler. There is no nested
% parallelism: one process worker executes one restart.
%
% The raw start exactly matches the corresponding frozen-SA restart start.
% Pattern Search then works directly in the linearly constrained physical
% constellation space. The final point is independently validated.

    taskTimer=tic;
    result=empty_result(identity,restartIndex,masterSeed,maxEvals);

    try
        if exist('patternsearch','file')~=2
            error('i13d_patternsearch_restart_task:ToolboxMissing', ...
                'Global Optimization Toolbox / patternsearch is unavailable on worker.');
        end

        Xraw=i13_shared_start_points(cfg,xBaseline,masterSeed,restartIndex);
        xRaw=Xraw(:,restartIndex);
        x0=AMI_functions.project_constellation_1D(xRaw,cfg);
        AMI_functions.assert_constellation_feasible(x0,cfg,1e-10);

        % Do NOT score x0 separately here. Pattern Search evaluates its
        % initial point inside MaxFunctionEvaluations. An extra logging
        % evaluation would give PS one more AMI call than the matched SA
        % budget. InitialFastAMI is therefore intentionally left NaN.
        initialFast=NaN;

        [A,b,Aeq,beq,lb,ub]=linear_constraints(cfg);

        opts=optimoptions('patternsearch', ...
            'Display',psPolicy.Display, ...
            'MaxFunctionEvaluations',double(maxEvals), ...
            'MaxIterations',double(maxEvals), ...
            'UseCompletePoll',logical(psPolicy.UseCompletePoll), ...
            'UseParallel',false, ...
            'MeshTolerance',double(psPolicy.MeshTolerance), ...
            'StepTolerance',double(psPolicy.StepTolerance), ...
            'FunctionTolerance',double(psPolicy.FunctionTolerance));

        fprintf('[I13D|%s] START | maxEvals=%d\n', ...
            char(identity.IdentityKey),maxEvals);

        optTimer=tic;
        [x,fval,exitflag,psOut]=patternsearch( ...
            @(z)-cfg.AMI_Evaluator(z(:)),x0,A,b,Aeq,beq,lb,ub,[],opts);
        optimizationSeconds=toc(optTimer);

        x=x(:);
        AMI_functions.assert_constellation_feasible(x,cfg,1e-8);
        fastMI=-fval;

        valTimer=tic;
        valMI=cfg.AMI_Validator(x);
        validationSeconds=toc(valTimer);

        result.success=true;
        result.x0Raw=xRaw;
        result.x0Feasible=x0;
        result.initialFastAMI=initialFast;
        result.xBest=x;
        result.bestMIFast=fastMI;
        result.bestMIValidated=valMI;
        result.validationGap=fastMI-valMI;
        result.optimizationSeconds=optimizationSeconds;
        result.validationSeconds=validationSeconds;
        result.taskWallSeconds=toc(taskTimer);
        result.exitflag=exitflag;
        if isfield(psOut,'funccount'),result.functionEvaluations=psOut.funccount;end
        if isfield(psOut,'iterations'),result.iterations=psOut.iterations;end
        if isfield(psOut,'meshsize'),result.finalMeshSize=psOut.meshsize;end
        result.output=psOut;

        checkpoint=struct();
        checkpoint.signatureText=signatureText;
        checkpoint.identity=identity;
        checkpoint.result=result;
        save_atomic(checkpointPath,checkpoint);

        fprintf('[I13D|%s] DONE | fast=%.6f | val=%.6f | delta=%+.3e | evals=%g | %.2fs\n', ...
            char(identity.IdentityKey),fastMI,valMI,fastMI-valMI, ...
            result.functionEvaluations,result.taskWallSeconds);

    catch ME
        result.success=false;
        result.taskWallSeconds=toc(taskTimer);
        result.errorIdentifier=string(ME.identifier);
        result.errorMessage=string(ME.message);
        result.errorReport=string(getReport(ME,'extended','hyperlinks','off'));
        fprintf(2,'[I13D|%s] FAILED | %s\n',char(identity.IdentityKey),ME.message);
    end
end

function r=empty_result(identity,restartIndex,masterSeed,maxEvals)
    r=struct();
    r.success=false;
    r.IdentityKey=string(identity.IdentityKey);
    r.M=identity.M;
    r.SNRdB=identity.SNRdB;
    r.SigmaR2=identity.SigmaR2;
    r.Replicate=identity.Replicate;
    r.Restart=restartIndex;
    r.masterSeed=double(masterSeed);
    r.maxEvals=maxEvals;
    r.x0Raw=[];
    r.x0Feasible=[];
    r.initialFastAMI=NaN;
    r.xBest=[];
    r.bestMIFast=NaN;
    r.bestMIValidated=NaN;
    r.validationGap=NaN;
    r.functionEvaluations=NaN;
    r.iterations=NaN;
    r.exitflag=NaN;
    r.finalMeshSize=NaN;
    r.optimizationSeconds=NaN;
    r.validationSeconds=NaN;
    r.taskWallSeconds=NaN;
    r.output=struct();
    r.errorIdentifier="";
    r.errorMessage="";
    r.errorReport="";
end

function [A,b,Aeq,beq,lb,ub]=linear_constraints(cfg)
    M=cfg.M;
    d=cfg.SA.minGap;

    A=zeros(M-1,M);
    for i=1:M-1
        A(i,i)=1;
        A(i,i+1)=-1;
    end
    b=-d*ones(M-1,1);

    nEq=1+double(cfg.SA.pinZero);
    Aeq=zeros(nEq,M);
    beq=zeros(nEq,1);
    Aeq(1,:)=1;
    beq(1)=M*cfg.P_avg;
    if cfg.SA.pinZero
        Aeq(2,1)=1;
        beq(2)=0;
    end

    lb=zeros(M,1);
    ub=[];
end

function save_atomic(path,checkpoint)
    folder=fileparts(path);
    if exist(folder,'dir')~=7,mkdir(folder);end
    tmp=[tempname(folder) '.mat'];
    cleanup=onCleanup(@() cleanup_tmp(tmp)); %#ok<NASGU>
    save(tmp,'checkpoint','-v7');
    [ok,msg]=movefile(tmp,path,'f');
    if ~ok,error('i13d_patternsearch_restart_task:CheckpointSave','%s',msg);end
end

function cleanup_tmp(path)
    if exist(path,'file')==2
        try,delete(path);catch,end
    end
end
