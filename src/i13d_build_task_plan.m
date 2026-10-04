function [casePlan,taskPlan] = i13d_build_task_plan(caseTable,F)
%I13D_BUILD_TASK_PLAN  Build the frozen 32-case / 224-restart I.13D plan.

    need={'M','SigmaX2','SNRdB','Replicate','Seed','NStarts', ...
        'SelectedDt','AMIValidated','GainBits'};
    missing=setdiff(need,caseTable.Properties.VariableNames);
    if ~isempty(missing)
        error('i13d_build_task_plan:Columns','Missing case-table columns: %s',strjoin(missing,', '));
    end

    keep=ismember(caseTable.M,F.MVec) & ...
        ismember(caseTable.SNRdB,F.SNRVec) & ...
        ismember(caseTable.Replicate,F.Replicates);
    sigKeep=false(height(caseTable),1);
    for s=F.TurbulenceVec
        sigKeep=sigKeep | abs(caseTable.SigmaX2-s)<1e-12;
    end
    T=caseTable(keep & sigKeep,:);
    T=sortrows(T,{'M','SNRdB','SigmaX2','Replicate'});

    if height(T)~=F.nOptimizationCases
        error('i13d_build_task_plan:CaseCount', ...
            'Expected %d cases, found %d.',F.nOptimizationCases,height(T));
    end

    CaseIndex=(1:height(T)).';
    SigmaR2=T.SigmaX2;
    PAMValidated=T.AMIValidated-T.GainBits;
    MaxEvalsPerStart=zeros(height(T),1);
    MaxIterSA=zeros(height(T),1);

    for c=1:height(T)
        iM=find(F.MVec==T.M(c),1);
        expectedStarts=F.RestartsByM(iM);
        if T.NStarts(c)~=expectedStarts
            error('i13d_build_task_plan:RestartCount', ...
                'Frozen SA restart count mismatch for M=%d.',T.M(c));
        end
        MaxEvalsPerStart(c)=F.MaxEvalsPerStart(iM);
        MaxIterSA(c)=F.MaxIterByM(iM);
    end

    casePlan=table(CaseIndex,T.M,T.SNRdB,SigmaR2,T.Replicate,T.Seed,T.NStarts, ...
        T.SelectedDt,MaxIterSA,MaxEvalsPerStart,PAMValidated,T.AMIValidated, ...
        'VariableNames',{'CaseIndex','M','SNRdB','SigmaR2','Replicate','Seed', ...
        'NStarts','SelectedDt','MaxIterSA','MaxEvalsPerStart', ...
        'PAMValidated','SAValidated'});

    TaskIndex=zeros(F.nRestartTasks,1);
    TaskCaseIndex=zeros(F.nRestartTasks,1);
    Restart=zeros(F.nRestartTasks,1);
    Priority=zeros(F.nRestartTasks,1);
    IdentityKey=strings(F.nRestartTasks,1);
    q=0;

    baseDt=0.01;
    for c=1:height(casePlan)
        R=casePlan(c,:);
        refined=R.SelectedDt < baseDt-1e-12;
        turbulent=R.SigmaR2 > 0;
        priority = 1e12*double(refined) + ...
                   1e10*double(turbulent) + ...
                   1e7*R.M + 1e4*R.SNRdB + 1e3*R.SigmaR2;
        for r=1:R.NStarts
            q=q+1;
            TaskIndex(q)=q;
            TaskCaseIndex(q)=c;
            Restart(q)=r;
            Priority(q)=priority;
            IdentityKey(q)=sprintf('M%d_SNR%d_sig%03d_rep%d_r%02d', ...
                R.M,round(R.SNRdB),round(1000*R.SigmaR2),R.Replicate,r);
        end
    end
    if q~=F.nRestartTasks
        error('i13d_build_task_plan:TaskCount','Expected %d tasks, built %d.',F.nRestartTasks,q);
    end

    taskPlan=table(TaskIndex,TaskCaseIndex,Restart,Priority,IdentityKey);
    taskPlan=sortrows(taskPlan,{'Priority','TaskIndex'},{'descend','ascend'});

    if numel(unique(taskPlan.IdentityKey))~=height(taskPlan)
        error('i13d_build_task_plan:DuplicateIdentity','Task identities are not unique.');
    end
end
