function S = commit_m_optimizer_statistics(T)
%COMMIT_M_OPTIMIZER_STATISTICS  Descriptive statistics of ALL 12 restarts.
% Figure 7 uses BestValidatedAMI (MAX), NEVER mean of replicate winners.
% Sample variance/std use denominator N-1. No uncertainty bars in Figure 7.
    assert(istable(T) && ...
        all(ismember({'Method','Restart','ValidatedAMI','FunctionEvaluations'}, ...
        T.Properties.VariableNames)), ...
        'CommitM:StatsSchema','Missing restart table columns.');
    methods=["GPS";"SA"];
    Method=methods;N=zeros(2,1);BestRestart=zeros(2,1);
    BestValidatedAMI=zeros(2,1);MeanValidatedAMI=zeros(2,1);
    VarianceValidatedAMI=zeros(2,1);StdValidatedAMI=zeros(2,1);
    MinValidatedAMI=zeros(2,1);MeanFunctionEvaluations=zeros(2,1);
    MinFunctionEvaluations=zeros(2,1);MaxFunctionEvaluations=zeros(2,1);
    for i=1:2
        R=sortrows(T(string(T.Method)==methods(i),:),'Restart');
        assert(height(R)==12 && isequal(double(R.Restart(:)),(1:12)') && ...
            all(isfinite(R.ValidatedAMI)) && all(isfinite(R.FunctionEvaluations)), ...
            'CommitM:StatsCount','Expected 12 finite restarts, 1:12, for %s.',methods(i));
        v=double(R.ValidatedAMI(:));
        [BestValidatedAMI(i),ix]=max(v);
        BestRestart(i)=R.Restart(ix);
        N(i)=numel(v);
        MeanValidatedAMI(i)=mean(v);
        VarianceValidatedAMI(i)=var(v,0);
        StdValidatedAMI(i)=std(v,0);
        MinValidatedAMI(i)=min(v);
        MeanFunctionEvaluations(i)=mean(R.FunctionEvaluations);
        MinFunctionEvaluations(i)=min(R.FunctionEvaluations);
        MaxFunctionEvaluations(i)=max(R.FunctionEvaluations);
    end
    S=table(Method,N,BestRestart,BestValidatedAMI, ...
        MeanValidatedAMI,VarianceValidatedAMI,StdValidatedAMI, ...
        MinValidatedAMI,MeanFunctionEvaluations, ...
        MinFunctionEvaluations,MaxFunctionEvaluations);
end
