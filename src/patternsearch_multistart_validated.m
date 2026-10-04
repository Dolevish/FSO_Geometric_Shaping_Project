function [out,runs] = patternsearch_multistart_validated(cfg,xBaseline,masterSeed,nStarts,maxEvalsPerStart,varargin)
%PATTERNSEARCH_MULTISTART_VALIDATED  Fair alternative-optimizer benchmark.
%
% Reviewer-2 utility. Runs MATLAB Pattern Search from the SAME raw restart
% starts used by SA, under the SAME feasible set and fast AMI objective.
% Every per-start winner is independently validated; the final winner is
% selected by validated AMI.
%
% The function requires Global Optimization Toolbox (patternsearch).
%
% Fairness controls:
%   - identical x0 generation to SA (i13_shared_start_points);
%   - same cfg.AMI_Evaluator and cfg.AMI_Validator;
%   - same non-negativity, x1=0, mean-power and d_min constraints;
%   - explicit objective-evaluation budget per start.
%
% maxEvalsPerStart counts Pattern Search function evaluations. For exact
% production-budget matching use SA maxIter+1 because each SA restart scores
% its initial state once and then one proposal per iteration.
%
% This function is a benchmark, not a claim of Pattern Search optimality.

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'Display','off',@text_scalar);
    addParameter(p,'UseCompletePoll',true,@logical_scalar);
    addParameter(p,'MeshTolerance',1e-6,@positive_scalar);
    addParameter(p,'StepTolerance',1e-8,@positive_scalar);
    addParameter(p,'FunctionTolerance',1e-8,@positive_scalar);
    parse(p,varargin{:}); o=p.Results;

    if exist('patternsearch','file')~=2
        error('patternsearch_multistart_validated:ToolboxMissing', ...
            ['patternsearch is unavailable. Reviewer-2 benchmark requires ' ...
             'MATLAB Global Optimization Toolbox.']);
    end
    if ~isfield(cfg,'AMI_Evaluator') || ~isfield(cfg,'AMI_Validator')
        error('patternsearch_multistart_validated:Evaluators','Both AMI evaluators are required.');
    end

    validateattributes(nStarts,{'numeric'},{'scalar','integer','positive','finite'});
    validateattributes(maxEvalsPerStart,{'numeric'},{'scalar','integer','>=',2,'finite'});

    Xraw=i13_shared_start_points(cfg,xBaseline,masterSeed,nStarts);
    [A,b,Aeq,beq,lb,ub]=linear_constraints(cfg);

    template=struct( ...
        'startIndex',NaN,'masterSeed',uint32(masterSeed), ...
        'x0_raw',[],'x0_feasible',[],'initialFastAMI',NaN, ...
        'x_best',[],'bestMIFast',NaN,'bestMIValidated',NaN, ...
        'validationGap',NaN,'runtime',NaN,'validationRuntime',NaN, ...
        'functionEvaluations',NaN,'iterations',NaN,'exitflag',NaN, ...
        'output',struct());
    runs=repmat(template,nStarts,1);

    for k=1:nStarts
        x0=AMI_functions.project_constellation_1D(Xraw(:,k),cfg);
        initialFast=cfg.AMI_Evaluator(x0);

        opts=optimoptions('patternsearch', ...
            'Display',char(string(o.Display)), ...
            'MaxFunctionEvaluations',double(maxEvalsPerStart), ...
            'MaxIterations',double(maxEvalsPerStart), ...
            'UseCompletePoll',logical(o.UseCompletePoll), ...
            'MeshTolerance',double(o.MeshTolerance), ...
            'StepTolerance',double(o.StepTolerance), ...
            'FunctionTolerance',double(o.FunctionTolerance));

        label=sprintf('[I13|PS]-%02d',k);
        fprintf('%s START | initMI=%.6f | maxEvals=%d\n',label,initialFast,maxEvalsPerStart);
        t=tic;
        [x,fval,exitflag,psOut]=patternsearch( ...
            @(z)-cfg.AMI_Evaluator(z(:)),x0,A,b,Aeq,beq,lb,ub,[],opts);
        runSeconds=toc(t);

        x=x(:);
        AMI_functions.assert_constellation_feasible(x,cfg,1e-8);
        fastMI=-fval;

        tv=tic;
        valMI=cfg.AMI_Validator(x);
        valSeconds=toc(tv);

        runs(k).startIndex=k;
        runs(k).masterSeed=uint32(masterSeed);
        runs(k).x0_raw=Xraw(:,k);
        runs(k).x0_feasible=x0;
        runs(k).initialFastAMI=initialFast;
        runs(k).x_best=x;
        runs(k).bestMIFast=fastMI;
        runs(k).bestMIValidated=valMI;
        runs(k).validationGap=fastMI-valMI;
        runs(k).runtime=runSeconds;
        runs(k).validationRuntime=valSeconds;
        if isfield(psOut,'funccount'),runs(k).functionEvaluations=psOut.funccount;end
        if isfield(psOut,'iterations'),runs(k).iterations=psOut.iterations;end
        runs(k).exitflag=exitflag;
        runs(k).output=psOut;

        fprintf('%s DONE | fast=%.6f | val=%.6f | delta=%+.3e | evals=%g | %.2fs\n', ...
            label,fastMI,valMI,fastMI-valMI,runs(k).functionEvaluations,runSeconds);
    end

    vals=[runs.bestMIValidated];
    fast=[runs.bestMIFast];
    [bestVal,idxVal]=max(vals);
    [bestFast,idxFast]=max(fast);

    out=struct();
    out.optimizer='MATLAB patternsearch';
    out.masterSeed=uint32(masterSeed);
    out.nStarts=nStarts;
    out.maxEvalsPerStart=maxEvalsPerStart;
    out.bestStartValidated=idxVal;
    out.bestMIValidated=bestVal;
    out.bestXValidated=runs(idxVal).x_best;
    out.bestStartFast=idxFast;
    out.bestMIFast=bestFast;
    out.selectionChanged=idxVal~=idxFast;
    out.totalFunctionEvaluations=sum([runs.functionEvaluations],'omitnan');
    out.totalRuntime=sum([runs.runtime]);
    out.totalValidationRuntime=sum([runs.validationRuntime]);
    out.maxAbsFastValidationGap=max(abs([runs.validationGap]));
    out.globalOptimumClaim=false;
end

function [A,b,Aeq,beq,lb,ub]=linear_constraints(cfg)
    M=cfg.M;
    d=cfg.SA.minGap;

    % x_i - x_{i+1} <= -d
    A=zeros(M-1,M);
    for i=1:M-1
        A(i,i)=1;
        A(i,i+1)=-1;
    end
    b=-d*ones(M-1,1);

    eqRows=1;
    if cfg.SA.pinZero,eqRows=2;end
    Aeq=zeros(eqRows,M); beq=zeros(eqRows,1);
    Aeq(1,:)=1;
    beq(1)=M*cfg.P_avg;
    if cfg.SA.pinZero
        Aeq(2,1)=1;
        beq(2)=0;
    end

    lb=zeros(M,1);
    ub=[];
end

function tf=text_scalar(v),tf=ischar(v)||(isstring(v)&&isscalar(v));end
function tf=logical_scalar(v),tf=islogical(v)&&isscalar(v);end
function tf=positive_scalar(v),tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>0;end
