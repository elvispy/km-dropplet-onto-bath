function theta = theta_max_radius(fn, fn_prime, fn_primeprime)
    % Calculates the theta that corresponds to the maximum radius possible
    % Using newton raphson with the equation fn'(theta) = 0
    % where fn(theta) = sin(theta) (R + zeta) [fn = zeta]
    
    
  
    % Derivative of the funtion
    fn_prime_anonymous = @(ang) sin(ang) .* fn_prime(ang) + cos(ang) .* (1+ fn(ang));
    
    % second derivative of the function
    fn_primeprime_anonymous = @(ang) 2*cos(ang) .* fn_prime(ang) + sin(ang) .* fn_primeprime(ang) - sin(ang) .* (1 + fn(ang));
    
    %We want a solution that satisfy theta > pi/2!
    theta = pi/2 + pi/4;
    tol_theta = 1e-10;
    n = 1;
    while abs(fn_prime_anonymous(theta)) >= tol_theta && n < 100
       theta = mod(theta - fn_prime_anonymous(theta)/fn_primeprime_anonymous(theta), pi);
       n = n + 1;
    end
    %disp(n);
    %assert(theta > pi/2-pi/6);

end