function [xFeasible,meta] = i13d_finalize_patternsearch_solution(xSolver,cfg,varargin)
%I13D_FINALIZE_PATTERNSEARCH_SOLUTION
% Canonicalize tiny numerical feasibility residuals in a Pattern Search result.
%
% Pattern Search enforces the linear constraints numerically and can return
% a point that misses an active spacing constraint by about solver tolerance
% (for example 0.00999903 when d_min=0.01).  The scientific feasible set of
% the project is defined by AMI_functions.project_constellation_1D, so I.13D
% applies that canonical projector once AFTER the optimizer terminates.
%
% This post-processing does not feed back into Pattern Search and therefore
% gives it no additional search evaluations.  A separate post-hoc fast AMI
% audit and the independent validator may then score xFeasible.
%
% The repair is accepted only when it is numerically tiny. A material repair
% fails loudly rather than silently changing the optimizer solution.

    p=inputParser; p.FunctionName=mfilename;
    addParameter(p,'MaxRepairInfNorm',1e-4,@positive_scalar);
    parse(p,varargin{:}); o=p.Results;

    xSolver=real(xSolver(:));
    if numel(xSolver)~=cfg.M || any(~isfinite(xSolver))
        error('i13d_finalize_patternsearch_solution:BadInput', ...
            'Pattern Search solution is non-finite or has wrong cardinality.');
    end

    before=constraint_diagnostics(xSolver,cfg);
    xFeasible=AMI_functions.project_constellation_1D(xSolver,cfg);
    repair=xFeasible-xSolver;

    meta=struct();
    meta.maxRepairInfNormAllowed=double(o.MaxRepairInfNorm);
    meta.repairInfNorm=norm(repair,inf);
    meta.repairL2Norm=norm(repair,2);
    meta.before=before;
    meta.after=constraint_diagnostics(xFeasible,cfg);
    meta.wasRepaired=meta.repairInfNorm>10*eps(max(1,norm(xSolver,inf)));

    if meta.repairInfNorm > o.MaxRepairInfNorm
        error('i13d_finalize_patternsearch_solution:MaterialRepair', ...
            ['Pattern Search output requires a material canonical repair: ' ...
             '||dx||_inf=%.3e > %.3e. Pre-repair min-gap deficit=%.3e, ' ...
             'mean error=%.3e, zero-anchor error=%.3e.'], ...
            meta.repairInfNorm,o.MaxRepairInfNorm,before.minGapDeficit, ...
            before.meanError,before.zeroAnchorError);
    end

    AMI_functions.assert_constellation_feasible(xFeasible,cfg,1e-10);
end


function d=constraint_diagnostics(x,cfg)
    gapReq=max(0,double(cfg.SA.minGap));
    if numel(x)>1
        minGapActual=min(diff(x));
    else
        minGapActual=Inf;
    end

    d=struct();
    d.minLevel=min(x);
    d.minGapActual=minGapActual;
    d.minGapRequired=gapReq;
    d.minGapDeficit=max(0,gapReq-minGapActual);
    d.meanActual=mean(x);
    d.meanTarget=cfg.P_avg;
    d.meanError=mean(x)-cfg.P_avg;
    if cfg.SA.pinZero
        d.zeroAnchorError=x(1);
    else
        d.zeroAnchorError=0;
    end
end

function tf=positive_scalar(v)
    tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&isreal(v)&&v>0;
end
