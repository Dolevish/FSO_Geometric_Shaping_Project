function [nStarts,meta] = i11_restart_policy(M)
%I11_RESTART_POLICY  Frozen restart-count policy for I.11 production.
%
% Production decision:
%   M = 4,8   -> 6 restarts
%   M >= 16   -> 8 restarts
%
% This is a conservative engineering freeze.  It is not a proof that the
% selected number is minimal and it is not a global-optimum guarantee.

    validateattributes(M,{'numeric'},{'scalar','integer','finite','>=',2},mfilename,'M');
    M=double(M);

    if M <= 8
        nStarts=6;
        reason='I11-v1 low-order production budget: M<=8 -> 6 restarts';
    else
        nStarts=8;
        reason='I11-v1 high-order conservative production budget: M>=16 -> 8 restarts';
    end

    meta=struct();
    meta.ruleVersion='I11-v1';
    meta.M=M;
    meta.nStarts=nStarts;
    meta.reason=reason;
    meta.note=['Measured/conservative production policy; not a proof of the ' ...
        'minimum restart count and not a global-optimum guarantee.'];
end
