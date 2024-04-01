function tann = calculate_tan(distance_to_axis, amplitudes)


    if isstruct(amplitudes); amplitudes = amplitudes.deformation_amplitudes; end
    angle = theta_from_cylindrical(distance_to_axis, amplitudes);
    zeta = zeta_generator(amplitudes);
    if size(amplitudes, 2) > 1; amplitudes = amplitudes'; end
    der = @(theta) sum(amplitudes .* collectdnPl(length(amplitudes), cos(theta)), 1);

    dzdr = @(theta)  (-sin(theta) .* (1 + zeta(theta)) - cos(theta) .* sin(theta) .* der(theta)) ./ ...
        (cos(theta) .* (1 + zeta(theta)) - sin(theta).^2 .* der(theta));
        
    tann =  tan(dzdr(angle));

end