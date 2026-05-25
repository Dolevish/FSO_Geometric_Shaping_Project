function [constellation_points, symbol_probs] = define_constellation_2D(M, P_avg, type)
% Defines different types of 2D complex constellations for phase noise channel.
%
% Inputs:
%   M     - Number of constellation points (e.g., 8)
%   P_avg - Desired average power of the constellation
%   type  - String: 'RANDOM_2D' (for SA start), 'PSK' (for testing)
%
% Outputs:
%   constellation_points - Row vector of COMPLEX constellation points [x1, x2, ..., xM]
%   symbol_probs         - Row vector of probabilities for each symbol [p1, p2, ..., pM]

    if nargin < 3
        type = 'RANDOM_2D'; % Default to random
    end

    % Assume equal probability for Geometric Shaping
    symbol_probs = (1/M) * ones(1, M);

    switch upper(type)
        case 'RANDOM_2D'
            temp_points = (randn(1, M) + 1i * randn(1, M));
          
            constellation_points = Enforce_Power_Constraint(temp_points, P_avg);

        case 'PSK'
            % Standard M-PSK constellation 
            angles = (0:M-1) * (2*pi/M);
            temp_points = exp(1i * angles);
            
            % scale PSK to P_avg
            constellation_points = temp_points * sqrt(P_avg);

        otherwise
            error('Unknown constellation type specified: %s', type);
    end

end