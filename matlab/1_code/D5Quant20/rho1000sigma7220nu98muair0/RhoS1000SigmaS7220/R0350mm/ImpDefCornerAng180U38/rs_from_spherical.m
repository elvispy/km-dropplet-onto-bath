function rs = rs_from_spherical(theta, amplitudes) % check pointwise function
    
    % This function returns the z-position of a point with respect to the
    % center of mass, with inputs in cylindrical coordinates.
    
    if isstruct(amplitudes); amplitudes = amplitudes.deformation_amplitudes; end
    if size(amplitudes, 2) > 1; amplitudes = amplitudes'; end
    
    zeta = zeta_generator(amplitudes);

    rs = sin(theta) .* (1 + zeta(theta)); % adimensionalized radius of the ball
   
end