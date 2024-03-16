function r = r_from_spherical(angles, amplitudes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% r_from_spherical.m - Returns cylindrical coordinates from spherical coordinates
%
% Calculate the radius in cylindrical coordinates that
% corresponds to the azimutal angle in spherical coordinates
% given certain amplitudes. (angle pi points downwards)
% r and angle are related by the relation
% r = sin(angle * (1 + \sum_{i=1}^{N}  A_l(i) * Pl(cos(angle)))
% where Pl is the l-th legendre Polynomial. (see
% https://en.wikipedia.org/wiki/Legendre_polynomials#Rodrigues'_formula_and_other_explicit_formulas) 

% Arguments:
% - angles: One dimensional scalar array such that size(r) = size(angles)
% angles(i) corresponds to the angle that gives radius r(i).
% - A_l: One dimensional scalar array or STRUCT with field
% "deformation_amplitudes"
%
% Outputs:
% - r: One dimensional scalar array of non-negative entries
%
% EXAMPLES
% r_from_spherical(pi, rand(10, 1));   % Returns 0.0
%
% Written by: Elvis Aguero- 01/01/2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isstruct(amplitudes); amplitudes = amplitudes.deformation_amplitudes; end
    if size(amplitudes, 2) > 1; amplitudes = amplitudes'; end
    r = sin(angles) .* (1 + sum(amplitudes .* collectPl(length(amplitudes), cos(angles)), 1));
%     r = arrayfun(@(angle) ...
%         sin(angle) * (1 + sum(dot(amplitudes, ...
%         arrayfun(@(idx) LP{idx}(cos(angle)), 1:d)))), theta);
end