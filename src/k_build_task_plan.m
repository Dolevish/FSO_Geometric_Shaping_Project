function [casePlan,taskPlan] = k_build_task_plan(K)
%K_BUILD_TASK_PLAN  Build Commit-K paired d_min cases and restart tasks.
%
% The case seed depends on (base seed, replicate, M, SNR, sigma_R^2) but NOT
% on d_min. Thus the five d_min values for a fixed (M,replicate) use the
% same master seed and restart substreams, giving a paired sensitivity study.

    nCases=K.nReplicateCases;
    CaseIndex=zeros(nCases,1);
    M=zeros(nCases,1);
    SNRdB=zeros(nCases,1);
    SigmaR2=zeros(nCases,1);
    MinGap=zeros(nCases,1);
    Replicate=zeros(nCases,1);
    Seed=zeros(nCases,1);
    NStarts=repmat(K.NStarts,nCases,1);
    MaxIter=repmat(K.MaxIter,nCases,1);
    SelectedDt=zeros(nCases,1);

    q=0;
    for m=K.MVec
        iM=find(K.MVec==m,1);
        if isfield(K,'MinGapByM')
            gapsForM=K.MinGapByM{iM};
        else
            gapsForM=K.MinGapVec;
        end
        for d=gapsForM
            for rep=K.Replicates
                q=q+1;
                repBaseSeed=double(K.baseSeed)+(rep-1)*10000019;
                seed=fso_result_utils.case_seed(repBaseSeed,m,K.SNRdB,K.SigmaR2);

                CaseIndex(q)=q;
                M(q)=m;
                SNRdB(q)=K.SNRdB;
                SigmaR2(q)=K.SigmaR2;
                MinGap(q)=d;
                Replicate(q)=rep;
                Seed(q)=seed;
                SelectedDt(q)=K.SelectedDtByM(iM);
            end
        end
    end
    if q~=nCases,error('k_build_task_plan:CaseCount','Internal case-count mismatch.');end

    casePlan=table(CaseIndex,M,SNRdB,SigmaR2,MinGap,Replicate,Seed, ...
        NStarts,MaxIter,SelectedDt);

    nTasks=K.nRestartTasks;
    TaskIndex=zeros(nTasks,1);
    TaskCaseIndex=zeros(nTasks,1);
    Restart=zeros(nTasks,1);
    Priority=zeros(nTasks,1);
    IdentityKey=strings(nTasks,1);
    q=0;

    maxGap=max(K.MinGapVec);
    for c=1:height(casePlan)
        R=casePlan(c,:);
        % M=32 first; within each M, smaller d_min first because it permits
        % the widest feasible amplitude range and can be computationally
        % more demanding under a candidate-adaptive y grid.
        priority=1e6*R.M + 1e4*(maxGap-R.MinGap);
        gapToken=strrep(sprintf('%.3f',R.MinGap),'.','p');
        for r=1:R.NStarts
            q=q+1;
            TaskIndex(q)=q;
            TaskCaseIndex(q)=c;
            Restart(q)=r;
            Priority(q)=priority;
            IdentityKey(q)=sprintf( ...
                'K_M%d_SNR%d_sig%03d_gap%s_rep%d_r%02d', ...
                R.M,round(R.SNRdB),round(1000*R.SigmaR2),gapToken,R.Replicate,r);
        end
    end
    if q~=nTasks,error('k_build_task_plan:TaskCount','Internal task-count mismatch.');end

    taskPlan=table(TaskIndex,TaskCaseIndex,Restart,Priority,IdentityKey);
    taskPlan=sortrows(taskPlan,{'Priority','TaskIndex'},{'descend','ascend'});

    if numel(unique(taskPlan.IdentityKey))~=height(taskPlan)
        error('k_build_task_plan:DuplicateIdentity','Task identities are not unique.');
    end
end
