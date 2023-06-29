function coefs = custom_project_amplitudes(angles, values, N, ~, ~) %, RelTol, AbsTol, PROBLEM_CONSTANTS, flag) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% custom_project_amplitudes.m - Projects function onto legendre polynomials
%
% Calculates the nth amplitude coefficient of a sampled function, where 
% N <= harmonics_qtt. Note that the zeroth mode is not calculated.
% THe following formula applies
% coefs(n) = (2n+1)/2 \int_{0}^{pi} ff(w) * sin(w) * P_l(cos(w)) dw
% where Pl is the l-th legendre Polynomial. (see
% https://en.wikipedia.org/wiki/Legendre_polynomials#Rodrigues'_formula_and_other_explicit_formulas)
% ff is the linear interpolation of the points given as arguments

% Arguments:
% - angles: angles at which the funciton is sampled
% - values: values of the sampled functions at angles given
% - N: number of harmonics to project the function to
% - PROBLEM_CONSTANTS: wether to artificially enhance linearity by adding
% points in between the angles
% - flag: Whether to activate enhancing
%
% Outputs:
% - coefs: Array of coefficients that correspond to the amplitudes.

% EXAMPLES
%
% Written by: Elvis Aguero- 27/06/2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Project from theta to x coordinate (in legendre terms)
cosang = cos(angles);

% Evaluate legendre polynomials once (from 0th order)
legendre_values = [ones(1, length(cosang)); collectPl(N+2, cosang)];

M = 3:2:(2*N+3); M = M';

% Calculate int P_n dx analytically (the privitive of P_n is proportional
% to P_{n+1} - P_{n-1}
Gn = (legendre_values(3:end, :) - legendre_values(1:(end-2), :)) ./M;
Gn = [cosang; Gn];
% Instantiate the integral int P_n dx from a to b
Gnab = diff(Gn, 1, 2);

% Calculate int x * P_n dx
intxPn = ((2:(N+1))'.* Gnab(3:end, :) + (1:N)' .* Gnab(1:N, :)) ./ M(1:N);

% Calculate a such that linear interpolation applies at the intervals
as = diff(values) ./ diff(cosang);
% Calculate independent term b such that a*x+b is the linear interpolator
bs = values(1:(end-1)) - as .* cosang(1:(end-1));

% Add up the contribution of every interval to the whole integral (negative
% because the orientation of the interval changes from theta to x
coefs = -M(1:(end-1))/2 .* sum(as .* intxPn + bs .* Gnab(2:(end-1), :), 2);

end