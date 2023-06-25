function rmax = maximum_contact_radius(amplitudes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% maximum_contact_radius.m - Calculates the maximum possible contact radius
%
% Given a dropplet described by the amplitudes, calculates the maximum
% possible contact radius, given by the first time the dropplet has a
% vertical tangent line. This maximum contact line is calculated using
% Newton-Raphson.
% The deformation is given by
% f(angle) = 1 +  \sum_{i=1}^{N}  amplitudes(i) * Pl(cos(angle)))
%
% Arguments:
% - amplitudes: One dimensional scalar array or STRUCT with field
% "deformation_amplitudes"
%
% Outputs:
% - rmax: Scalar value. 
%
% EXAMPLES
% maximum_contact_radius([0, 0]);   % Returns 1
%
% Written by: Elvis Aguero- 01/01/2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if isstruct(amplitudes); amplitudes = amplitudes.deformation_amplitudes; end
    if size(amplitudes, 2) > 1; amplitudes = amplitudes'; end
    %collectdnPl = PROBLEM_CONSTANTS.collectdnPl;
    %collectd2nPl = PROBLEM_CONSTANTS.collectd2nPl;
    order = length(amplitudes);
    zeta = zeta_generator(amplitudes); 
    drdtheta = @(theta) cos(theta) * (1 + zeta(theta)) - ...
        sin(theta)^2 * sum(times(amplitudes, collectdnPl(order, cos(theta))));

    dr2dtheta2 = @(theta) - sin(theta) * (1 + zeta(theta)) - 2 * cos(theta) * sin(theta) * ...
        sum(times(amplitudes, collectdnPl(order, cos(theta))), 1) + ...
        sin(theta)^3 .* sum(times(amplitudes, collectdnPl(order, cos(theta), 2)), 1);

    theta = pi/2 + pi/4;
    tol_theta = 1e-7;
    n = 1;
    % Newton Method!
    while abs(drdtheta(theta)) >= tol_theta && n < 150
        theta = mod(theta - drdtheta(theta)/dr2dtheta2(theta) - 1e-4, pi) + 1e-4; % If solution is close to pi, theta is unstable with mod function (therefore 1e-4 added)
        n = n + 1;
        if n == 50
            theta = 3.14159/2;
        elseif n == 100
            theta = rand()/100 + pi/2;
        end
        if n== 149
            % We will calculate it manually
            theta_min = pi/3;
            theta_max = 4*pi/5;
            thetas = linspace(theta_min, theta_max, 1000);
            rmax = max(r_from_spherical(thetas, amplitudes));
            
            return
        end
    end

    rmax = r_from_spherical(theta, amplitudes);

end