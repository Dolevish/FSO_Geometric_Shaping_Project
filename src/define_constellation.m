function x = define_constellation(M, P_avg, imdd_mode, powerConstraint)
%DEFINE_CONSTELLATION  Build a normalized uniform PAM constellation.
%
%   x = define_constellation(M, P_avg, imdd_mode, powerConstraint)
%
% Inputs:
%   M               - constellation cardinality, integer >= 2
%   P_avg           - target power/intensity constraint, positive scalar
%   imdd_mode       - logical. true -> non-negative IM/DD PAM [0,...,M-1]
%                              false -> bipolar PAM [-(M-1),...,M-1]
%   powerConstraint - "mean" for average optical intensity,
%                     "meansquare" for average electrical power
%
% The IEEE FSO work uses:
%       define_constellation(M, P_avg, true, "mean")
%
% Strict input checking is intentional: passing "mean" as the third
% argument used to be silently interpreted as imdd_mode and caused the
% baseline to be normalized with the wrong power definition.

    if nargin < 3 || isempty(imdd_mode),       imdd_mode = false; end
    if nargin < 4 || isempty(powerConstraint), powerConstraint = "meansquare"; end

    validateattributes(M, {'numeric'}, ...
        {'scalar','integer','>=',2,'finite'}, mfilename, 'M', 1);
    validateattributes(P_avg, {'numeric'}, ...
        {'scalar','real','positive','finite'}, mfilename, 'P_avg', 2);

    if ~(islogical(imdd_mode) && isscalar(imdd_mode))
        error('define_constellation:BadIMDDMode', ...
            ['imdd_mode must be a logical scalar (true/false). ' ...
             'For this project use define_constellation(M,P_avg,true,"mean").']);
    end

    powerConstraint = string(powerConstraint);
    if ~isscalar(powerConstraint) || ~any(powerConstraint == ["mean","meansquare"])
        error('define_constellation:BadPowerConstraint', ...
            'powerConstraint must be "mean" or "meansquare".');
    end

    if imdd_mode
        % IM/DD intensity levels (non-negative, includes zero level).
        x = 0:(M-1);
    else
        % Bipolar M-PAM.
        x = -(M-1):2:(M-1);
    end

    x = Enforce_Power_Constraint(x(:), P_avg, powerConstraint);
    x = sort(x(:)).';
end
