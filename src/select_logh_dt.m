function [dt,meta] = select_logh_dt(M,SNRdB,sigmaX2,varargin)
%SELECT_LOGH_DT  Deterministic per-channel-case log-h spacing policy.
%
% Commit I.9 policy derived from the I.7/I.8 numerical diagnostics.
% The selected spacing depends ONLY on the physical channel case
% (M,SNR_dB,sigma_X^2). It never depends on the current constellation, SA
% iteration, restart, or random seed; therefore dt is constant throughout a
% complete SA optimization case.
%
% Default case-adaptive policy:
%   base dt    = 0.01
%   refined dt = 0.005
%
% Refinement rule:
%   SNR >= 30 dB AND
%      [ (M >= 32 AND sigma_X^2 >= 0.1) OR
%        (M >= 16 AND sigma_X^2 >= 0.3) ]
%
% On the current production grid M={4,8,16,32}, SNR<=30 dB and
% sigma_X^2={0,0.1,0.2,0.3}, this selects exactly the four I.8 stress cases:
%   (32,30,0.1), (32,30,0.2), (32,30,0.3), (16,30,0.3).
%
% The >= formulation is intentionally conservative for any future extension
% beyond the current grid: conditions at least as severe are refined rather
% than silently reverting to the coarse spacing.
%
% Name-value options:
%   'Policy'     'case-adaptive' (default) or 'fixed'
%   'BaseDt'     coarse spacing for case-adaptive mode (default 0.01)
%   'RefinedDt'  refined spacing for hard cases (default 0.005)
%   'FixedDt'    spacing when Policy='fixed' (default 0.01)
%
% Outputs:
%   dt    selected numerical spacing
%   meta  explicit reproducibility metadata including rule and reason

    p = inputParser;
    p.FunctionName = mfilename;
    addRequired(p,'M',@positive_integer);
    addRequired(p,'SNRdB',@finite_scalar);
    addRequired(p,'sigmaX2',@nonnegative_scalar);
    addParameter(p,'Policy','case-adaptive',@(v) ischar(v) || (isstring(v)&&isscalar(v)));
    addParameter(p,'BaseDt',0.01,@positive_scalar);
    addParameter(p,'RefinedDt',0.005,@positive_scalar);
    addParameter(p,'FixedDt',0.01,@positive_scalar);
    parse(p,M,SNRdB,sigmaX2,varargin{:});
    o = p.Results;

    policy = lower(string(o.Policy));
    if ~ismember(policy,["case-adaptive","fixed"])
        error('select_logh_dt:BadPolicy', ...
            'Policy must be ''case-adaptive'' or ''fixed'', got ''%s''.',policy);
    end
    if o.RefinedDt > o.BaseDt
        error('select_logh_dt:BadSpacingOrder', ...
            'RefinedDt (%.6g) must not exceed BaseDt (%.6g).',o.RefinedDt,o.BaseDt);
    end

    tol = 1e-12;
    hardCase = double(SNRdB) >= 30-tol && ...
        ((double(M) >= 32 && double(sigmaX2) >= 0.1-tol) || ...
         (double(M) >= 16 && double(sigmaX2) >= 0.3-tol));

    if policy == "fixed"
        dt = double(o.FixedDt);
        isRefined = false;
        reason = sprintf('fixed override: dt=%.6g',dt);
    elseif hardCase
        dt = double(o.RefinedDt);
        isRefined = true;
        reason = sprintf(['I9 high-SNR/high-order refinement: SNR>=30 and ' ...
            '((M>=32,sigma>=0.1) or (M>=16,sigma>=0.3)); dt=%.6g'],dt);
    else
        dt = double(o.BaseDt);
        isRefined = false;
        reason = sprintf('I9 base case: refinement rule not triggered; dt=%.6g',dt);
    end

    meta = struct();
    meta.policy = char(policy);
    meta.ruleVersion = 'I9-v1';
    meta.selectedDt = dt;
    meta.baseDt = double(o.BaseDt);
    meta.refinedDt = double(o.RefinedDt);
    meta.fixedDt = double(o.FixedDt);
    meta.isRefined = logical(isRefined);
    meta.hardCaseRuleTriggered = logical(hardCase);
    meta.M = double(M);
    meta.SNRdB = double(SNRdB);
    meta.sigmaX2 = double(sigmaX2);
    meta.reason = reason;
end

function tf=positive_scalar(v), tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&isreal(v)&&v>0; end
function tf=nonnegative_scalar(v), tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&isreal(v)&&v>=0; end
function tf=finite_scalar(v), tf=isnumeric(v)&&isscalar(v)&&isfinite(v)&&isreal(v); end
function tf=positive_integer(v), tf=positive_scalar(v)&&mod(v,1)==0; end
