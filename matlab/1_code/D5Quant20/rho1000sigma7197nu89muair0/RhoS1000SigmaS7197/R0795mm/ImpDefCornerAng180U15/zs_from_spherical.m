function zs = zs_from_spherical(theta, zeta) % check pointwise function
    
    % This function returns the z-position of a point with respect to the
    % center of mass, with inputs in cylindrical coordinates.
    if size(theta, 1) > 1; theta = theta'; end
    zs = cos(theta) .* (1 + zeta(theta)); % adimensionalized radius of the ball
    % s = cos(theta);
    
    % for ii = 1:length(amplitudes)
    %     zs = zs + amplitudes(ii) * legendreP(ii, s);
    % end
    % zs = zs .* s;
end