function S=commit_o_optimizer_statistics(T)
%COMMIT_O_OPTIMIZER_STATISTICS  All-12 descriptive results for 3 independent arms.
% Sample variance/std (N-1). MAX only is plotted for SA and GPS tight.
    required={'Method','Restart','ValidatedAMI','FunctionEvaluations'};
    assert(istable(T)&&all(ismember(required,T.Properties.VariableNames)), ...
        'CommitO:StatsColumns','Missing restart table columns.');
    Method=["GPS";"SA";"GPS_CONTROL"];
    n=numel(Method);
    N=zeros(n,1);BestRestart=zeros(n,1);BestValidatedAMI=zeros(n,1);
    MeanValidatedAMI=zeros(n,1);VarianceValidatedAMI=zeros(n,1);
    StdValidatedAMI=zeros(n,1);MinValidatedAMI=zeros(n,1);
    MeanFunctionEvaluations=zeros(n,1);MinFunctionEvaluations=zeros(n,1);
    MaxFunctionEvaluations=zeros(n,1);
    for i=1:n
        A=sortrows(T(string(T.Method)==Method(i),:),'Restart');
        assert(height(A)==12 && isequal(double(A.Restart(:)),(1:12)') && ...
            all(isfinite(A.ValidatedAMI)) && ...
            all(isfinite(A.FunctionEvaluations)), ...
            'CommitO:StatsCount','Expected 12 finite restarts for %s.',Method(i));
        x=double(A.ValidatedAMI(:));
        [BestValidatedAMI(i),ix]=max(x); BestRestart(i)=A.Restart(ix);
        N(i)=numel(x);MeanValidatedAMI(i)=mean(x);
        VarianceValidatedAMI(i)=var(x,0);StdValidatedAMI(i)=std(x,0);
        MinValidatedAMI(i)=min(x);
        MeanFunctionEvaluations(i)=mean(A.FunctionEvaluations);
        MinFunctionEvaluations(i)=min(A.FunctionEvaluations);
        MaxFunctionEvaluations(i)=max(A.FunctionEvaluations);
    end
    S=table(Method,N,BestRestart,BestValidatedAMI, ...
        MeanValidatedAMI,VarianceValidatedAMI,StdValidatedAMI, ...
        MinValidatedAMI,MeanFunctionEvaluations,MinFunctionEvaluations, ...
        MaxFunctionEvaluations);
end
