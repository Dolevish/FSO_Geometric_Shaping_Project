function [out,runs] = fminsearch_multistart_validated(cfg,xBaseline,masterSeed,nStarts,maxEvalsPerStart,varargin)
%FMINSEARCH_MULTISTART_VALIDATED  Nelder-Mead alternative optimizer benchmark.
%
% Reviewer-2 utility for MATLAB installations without Global Optimization
% Toolbox. Uses core-MATLAB fminsearch (Nelder-Mead) with the SAME canonical
% constellation projector used by SA:
%
%       objective(z) = -AMI_fast(Project(z)).
%
% Thus every AMI evaluation is performed on the same feasible set:
% nonnegative sorted levels, optional x1=0, d_min spacing and fixed mean
% optical power. The raw restart starts are exactly those used by SA.
%
% Every per-start winner is independently validated and final selection is
% by validated AMI. maxEvalsPerStart is an objective-evaluation budget.
%
% This is an optimizer-comparison tool, not a global-optimum claim.

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'Display','off',@text_scalar);
    addParameter(p,'TolX',1e-7,@positive_scalar);
    addParameter(p,'TolFun',1e-8,@positive_scalar);
    parse(p,varargin{:}); o=p.Results;

    if exist('fminsearch','file')~=2
        error('fminsearch_multistart_validated:Unavailable','fminsearch is unavailable.');
    end
    if ~isfield(cfg,'AMI_Evaluator') || ~isfield(cfg,'AMI_Validator')
        error('fminsearch_multistart_validated:Evaluators','Both AMI evaluators are required.');
    end

    validateattributes(nStarts,{'numeric'},{'scalar','integer','positive','finite'});
    validateattributes(maxEvalsPerStart,{'numeric'},{'scalar','integer','>=',2,'finite'});

    Xraw=i13_shared_start_points(cfg,xBaseline,masterSeed,nStarts);

    template=struct( ...
        'startIndex',NaN,'masterSeed',uint32(masterSeed), ...
        'x0_raw',[],'x0_feasible',[],'initialFastAMI',NaN, ...
        'searchVectorFinal',[],'x_best',[],'bestMIFast',NaN, ...
        'bestMIValidated',NaN,'validationGap',NaN, ...
        'runtime',NaN,'validationRuntime',NaN, ...
        'functionEvaluations',NaN,'iterations',NaN,'exitflag',NaN, ...
        'output',struct());
    runs=repmat(template,nStarts,1);

    opts=optimset('Display',char(string(o.Display)), ...
        'MaxFunEvals',double(maxEvalsPerStart), ...
        'MaxIter',double(maxEvalsPerStart), ...
        'TolX',double(o.TolX),'TolFun',double(o.TolFun));

    for k=1:nStarts
        x0=AMI_functions.project_constellation_1D(Xraw(:,k),cfg);
        initialFast=cfg.AMI_Evaluator(x0);

        label=sprintf('[I13|NM]-%02d',k);
        fprintf('%s START | initMI=%.6f | maxEvals=%d\n',label,initialFast,maxEvalsPerStart);
        t=tic;
        [z,fval,exitflag,nmOut]=fminsearch( ...
            @(u) projected_negative_ami(u,cfg),x0,opts);
        runSeconds=toc(t);

        x=AMI_functions.project_constellation_1D(z,cfg);
        AMI_functions.assert_constellation_feasible(x,cfg,1e-9);

        % Re-evaluate the final projected constellation explicitly so the
        % stored fast value corresponds exactly to the stored x_best.
        fastMI=cfg.AMI_Evaluator(x);
        tv=tic;
        valMI=cfg.AMI_Validator(x);
        valSeconds=toc(tv);

        runs(k).startIndex=k;
        runs(k).masterSeed=uint32(masterSeed);
        runs(k).x0_raw=Xraw(:,k);
        runs(k).x0_feasible=x0;
        runs(k).initialFastAMI=initialFast;
        runs(k).searchVectorFinal=z(:);
        runs(k).x_best=x(:);
        runs(k).bestMIFast=fastMI;
        runs(k).bestMIValidated=valMI;
        runs(k).validationGap=fastMI-valMI;
        runs(k).runtime=runSeconds;
        runs(k).validationRuntime=valSeconds;
        if isfield(nmOut,'funcCount'),runs(k).functionEvaluations=nmOut.funcCount;end
        if isfield(nmOut,'iterations'),runs(k).iterations=nmOut.iterations;end
        runs(k).exitflag=exitflag;
        runs(k).output=nmOut;

        fprintf('%s DONE | fast=%.6f | val=%.6f | delta=%+.3e | evals=%g | %.2fs\n', ...
            label,fastMI,valMI,fastMI-valMI,runs(k).functionEvaluations,runSeconds);
    end

    vals=[runs.bestMIValidated];
    fast=[runs.bestMIFast];
    [bestVal,idxVal]=max(vals);
    [bestFast,idxFast]=max(fast);

    out=struct();
    out.optimizer='MATLAB fminsearch (Nelder-Mead + canonical projection)';
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

function f=projected_negative_ami(z,cfg)
    x=AMI_functions.project_constellation_1D(z(:),cfg);
    f=-cfg.AMI_Evaluator(x);
    if ~(isscalar(f)&&isfinite(f))
        error('fminsearch_multistart_validated:InvalidObjective','AMI objective returned invalid value.');
    end
end

function tf=text_scalar(v),tf=ischar(v)||(isstring(v)&&isscalar(v));end
function tf=positive_scalar(v),tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>0;end
