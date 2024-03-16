function angle = theta_from_cylindrical(r, A_l)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% theta_from_cylindrical.m - Returns spherical coordinates from cylindrical coordinates
%
%  Calculate the azimutal angle in spherical coordinates
% that corresponds to radius r in cylindrical coordinates, 
% given amplitudes A_l. (angle pi points downwards)
% This algorithm uses Newton-Raphson to approximate the angle
% r and angle are related by the relation
% r = sin(angle) * (1 + \sum_{i=1}^{N}  A_l(i) * Pl(cos(angle)))
% where Pl is the l-th legendre Polynomial. (see
% https://en.wikipedia.org/wiki/Legendre_polynomials#Rodrigues'_formula_and_other_explicit_formulas) 

% Arguments:
% - r: One dimensional scalar array of non-negative entries
% - A_l: One dimensional scalar matrix or STRUCT with field
% "deformation_amplitudes"
%
% Outputs:
% - angles: One dimensional scalar array such that size(r) = size(angles)
% angles(i) corresponds to the angle that gives radius r(i).

% EXAMPLES
% theta_from_cylindrical(0, rand(10, 1));   % Returns pi
%
% Written by: Elvis Aguero- 01/01/2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if isstruct(A_l); A_l = A_l.deformation_amplitudes; end
    if size(A_l, 2) > 1; A_l = A_l'; end

    zeta = zeta_generator(A_l);
    
    % Derivative of the function
    f_prime = @(theta) cos(theta) .* (1 + zeta(theta)) - sin(theta).^2 .* sum(times(A_l, collectdnPl(length(A_l), cos(theta))), 1);
    
    
    angle = zeros(size(r));

    % Newton Method!
    for ii = 1:length(r)
        % Function to be minimized00
        f_objective = @(theta) sin(theta) .* (1 + zeta(theta)) - r(ii);

        theta = pi - 0.1;
        tol_theta = 1e-7;
        n = 1;
        
        while abs(f_objective(theta)) >= tol_theta && n < 150
            theta = mod(theta - f_objective(theta)/f_prime(theta) - 1e-4, pi/2) + 1e-4 + pi/2; % If solution is close to pi, theta is unstable with mod function (therefore 1e-4 added)
            n = n + 1;
            if n == 50
                theta = 3.14159;
            elseif n == 100
                theta = rand() * pi/2 + pi/2;
            end
        end
        angle(ii) = theta;
    end

end