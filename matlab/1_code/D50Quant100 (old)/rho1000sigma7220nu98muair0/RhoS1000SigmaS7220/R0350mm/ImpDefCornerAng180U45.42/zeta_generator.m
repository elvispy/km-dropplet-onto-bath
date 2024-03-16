function f = zeta_generator(amplitudes, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% zeta_generator.m - function handle of spectral coordinates given
% amplitudes
%
% Returns a function handle f, such that 
% f(angle) = \sum_{i=1}^{N}  amplitudes(i) * Pl(cos(angle)))
% where Pl is the l-th legendre Polynomial. (see
% https://en.wikipedia.org/wiki/Legendre_polynomials#Rodrigues'_formula_and_other_explicit_formulas) 
% The angle must satisfy 0 <= angle <= pi

% Arguments:
% - amplitudes: One dimensional scalar array or STRUCT with field
% "deformation_amplitudes"
% - varargin: Use of varargin is deprecated as of 01/23/2023
%
% Outputs:
% - f: function handle

% EXAMPLES
% zeta_generator([0 1]);   % Returns 1/2 * (3 * cos(angle).^2 - 1)
%
% Written by: Elvis Aguero- 01/01/2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if isstruct(amplitudes); amplitudes = amplitudes.deformation_amplitudes; end
    
    order = length(amplitudes);
    
    if nargin == 1
        if size(amplitudes, 2) > 1; amplitudes = amplitudes'; end
        f = @(theta) sum(amplitudes .* collectPl(order, cos(theta)), 1);
    else
        warning("I'm deprecated!");
        LP = varargin{1}.LEGENDRE_POLYNOMIALS;
        f = @(theta) arrayfun(@(ang) sum(times(amplitudes, arrayfun(@(idx) LP{idx}(cos(ang)), 1:order))), theta);
    end
end