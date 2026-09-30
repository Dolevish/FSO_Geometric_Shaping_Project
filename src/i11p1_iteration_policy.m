function [maxIter,meta] = i11p1_iteration_policy(M)
%I11P1_ITERATION_POLICY  Frozen per-M SA iteration budget for production.
%
% I.11.1 production policy:
%   M = 32      -> 8000 iterations/restart
%   all others  -> 4000 iterations/restart
%
% The larger M=32 budget is a conservative production choice. It is not a
% claim that 8000 iterations are minimal, optimal, or required in every
% physical channel condition.

    validateattributes(M,{'numeric'},{'scalar','integer','finite','>=',2},mfilename,'M');
    M=double(M);

    if M==32
        maxIter=8000;
        reason='I.11.1 conservative M=32 production budget';
    else
        maxIter=4000;
        reason='I.11.1 standard production budget for M ~= 32';
    end

    meta=struct('ruleVersion','I11.1-v1','M',M,'maxIter',maxIter, ...
        'reason',reason,'globalOptimumClaim',false);
end
