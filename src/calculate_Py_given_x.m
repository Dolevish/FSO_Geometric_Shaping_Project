function p_out = calculate_Py_given_x(y_in, x, params)
% Calculates the conditional probability density function (PDF) P(y|x)
% for the FSO turbulent channel with additive white Gaussian noise (AWGN).
% It uses a robust numerical integration method based on a change of variables.
%
% Mathematically, it calculates:
% P(y|x) = ∫[from h=0 to inf] P(y|x,h) * P(h) dh
% where:
%   P(y|x,h) = N(y; mean=R*x*h, variance=sigma_n^2)  (Gaussian noise PDF)
%   P(h)      = Log-normal PDF for weak turbulence, parameterized by sigma_X^2
%
% To improve numerical stability, the integration variable is changed to t = ln(h).
% Since h is log-normally distributed, t is normally distributed:
%   t ~ N(mu_t, sig_t^2)
%   where mu_t = -2*sigma_X^2 and sig_t^2 = 4*sigma_X^2
%
% The integral becomes:
% P(y|x) = ∫[from t=-inf to inf] N(y; R*x*e^t, sigma_n^2) * N(t; mu_t, sig_t^2) dt
%
% This function supports vector input `y_in` for efficiency when called by
% functions like 'integral' or 'quadgk' using 'ArrayValued',true.

    % --- Log-Normal Parameters for t = ln(h) ---
    % Correct parametrization: sigma_X_sq = Var[h]/E[h]^2 = scintillation index.
    % For h ~ LogNormal(mu_t, sig_t^2), the scintillation index gives:
    %   sig_t^2 = log(1 + sigma_X_sq)
    %   mu_t    = -0.5 * sig_t^2          (ensures E[h] = 1)
    % This is consistent with config_params.m and AMI_noCSI_GH.m.
    % NOTE: The old parametrization (mu_t=-2*sX2, sig_t=2*sqrt(sX2)) was INCORRECT
    %       and caused a large gap between fast evaluator and validation.
    sig_t_sq = log(1 + params.sigma_X_sq); % Variance of ln(h)
    sig_t    = sqrt(sig_t_sq);             % Standard deviation of ln(h)
    mu_t     = -0.5 * sig_t_sq;           % Mean of ln(h)  [ensures E[h]=1]

    % --- Handle AWGN Case Separately ---
    if sig_t < 1e-6
        % No turbulence: h = exp(mu_t) ≈ exp(0) = 1 (since mu_t → 0 as sigma_X_sq → 0)
        signal_mean    = params.R * x * exp(mu_t);   % ≈ R*x when turbulence→0
        sigma_n_sq_eff = params.sigma_n_sq;
        % Calculate the Gaussian PDF directly. Note the element-wise operations (./ .^ )
        % to handle vector y_in.
        p_out = (1./sqrt(2*pi*sigma_n_sq_eff)) .* exp( -(y_in - signal_mean).^2 ./ (2*sigma_n_sq_eff) );
        % Ensure the result is non-negative and finite (robustness check)
        p_out(p_out < 0 | ~isfinite(p_out)) = 0;
        return; % Exit the function early as calculation is done
    end

    % --- Define PDFs for Turbulent Case Integration ---

    % Define P(y|x,t) = N(y; R*x*e^t, sigma_n^2) as an anonymous function @(t)
    % This represents the probability of receiving y given transmitted x and a specific
    % turbulence state represented by t = ln(h).
    % 't' is scalar input, 'y_in' can be a vector. Element-wise ops ensure output matches y_in size.
    gauss_y_given_t = @(t) (1./sqrt(2*pi*params.sigma_n_sq)) .* ...
                           exp( -(y_in - params.R.*x.*exp(t)).^2 ./ (2*params.sigma_n_sq) );

    % Define P(t) = N(t; mu_t, sig_t^2) as an anonymous function @(t)
    % This is the PDF of the normally distributed variable t = ln(h).
    % 't' is scalar input, output is scalar.
    norm_t = @(t) (1/(sqrt(2*pi)*sig_t)) * exp( -(t-mu_t).^2./(2*sig_t_sq) );

    % Define the full integrand: P(y|x,t) * P(t)
    % This anonymous function @(t) calculates the product of the two PDFs defined above.
    % Since gauss_y_given_t returns a vector (size of y_in) and norm_t returns a scalar,
    % the result of the multiplication is a vector of the same size as y_in.
    % This is required for 'ArrayValued',true in the 'integral' function.
    integrand = @(t) gauss_y_given_t(t) .* norm_t(t);

    % --- Integration ---
    % Integrate the 'integrand' function with respect to 't'.
    % Instead of integrating from -inf to inf, which can be numerically unstable,
    % integrate over a finite range where most of the probability mass of P(t) lies.
    % A range of mu_t +/- 8*sig_t contains virtually all the probability for a normal distribution.
    integration_limit_stds = 8; % Number of standard deviations from the mean for integration limits
    lo = mu_t - integration_limit_stds * sig_t; % Lower integration limit
    hi = mu_t + integration_limit_stds * sig_t; % Upper integration limit

    try
        % Perform the numerical integration using MATLAB's 'integral' function.
        % 'integrand': the function to integrate.
        % 'lo', 'hi': the integration limits.
        % 'RelTol', 'AbsTol': specify the desired relative and absolute error tolerances.
        % 'ArrayValued', true: indicates that 'integrand' accepts scalar 't' and returns
        %                      an array (matching the size of y_in), enabling efficient
        %                      calculation for vector y_in.
        p = integral(integrand, lo, hi, 'RelTol',1e-9, 'AbsTol',1e-12, 'ArrayValued',true);
    catch ME
        % If the integration fails (e.g., doesn't converge, encounters NaN/Inf),
        % catch the error, display a warning, and return zero probability.
        warning('Integration failed in calculate_Py_given_x for x=%.2f. Error: %s', x, ME.message);
        % Return zeros matching the size of y_in if integration fails
        p = zeros(size(y_in));
    end

    % --- Cleanup and Output ---
    % Ensure the final probability values are non-negative and finite.
    % Numerical inaccuracies might sometimes produce very small negative numbers.
    p(p < 0 | ~isfinite(p)) = 0;
    % Assign the calculated probability vector to the output variable.
    p_out = p;

end % End of function calculate_Py_given_x