function coefs = custom_project_amplitudes(angles, values, N, RelTol, AbsTol, PROBLEM_CONSTANTS, flag) 
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
% - harmonics_qtt: number of harmonics to project the function to
% - RelTol: Relative Tolerance to be achieved
% - Abstol: Asolute Tolerance to be achieved
% - PROBLEM_CONSTANTS: Struct that has a field with the precomputed weights
% for the Curtis-Clensahw Quadrature.
% - flag: Boolean that flags whether to use Curtis-CLenshaw. If set to
% false, it uses Curtis Clenshaw.
%
% Outputs:
% - coefs: Array of coefficients that correspond to the amplitudes.

% EXAMPLES
%
% Written by: Elvis Aguero- 27/06/2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




end