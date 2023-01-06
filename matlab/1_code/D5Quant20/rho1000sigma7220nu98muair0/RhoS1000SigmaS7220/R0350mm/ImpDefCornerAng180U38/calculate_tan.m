function tann = calculate_tan(distance_to_axis, fn, fn_prime)
    % Calculates the tangent of the sphere at distance distance_to_axis(ii)
    
    
    M = length(distance_to_axis);
    thetas = zeros(M, 1);
    for ii =1:(M-2)
        if ii > 1
            thetas(ii) = theta_from_cylindrical(distance_to_axis(ii), fn, fn_prime, thetas(ii-1));
        else
            thetas(ii) = theta_from_cylindrical(distance_to_axis(ii), fn, fn_prime, pi);
        end
    end
    RR = arrayfun(fn, thetas);
    RRprime = arrayfun(fn_prime, thetas); 
    
    tann = double((-sin(thetas) .* (1 + RR) + cos(thetas).* RRprime) ./ ...
            (cos(thetas) .* (1 + RR) + sin(thetas) .* RRprime));
end