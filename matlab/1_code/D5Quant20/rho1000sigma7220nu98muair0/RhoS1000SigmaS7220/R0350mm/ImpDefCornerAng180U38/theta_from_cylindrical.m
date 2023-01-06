
function theta = theta_from_cylindrical(r, fn, fn_prime, guess) 
    % Calculates the theta that corresponds to a radius by Newton-Raphson
    % Where fn(theta) = sin(theta) * (R + zeta)
    %N = length(amplitudes); syms x;
    
    %if size(amplitudes, 1) > 1; amplitudes = amplitudes'; end
    % Derivative of the function i'm trying to find the zero of
    fn_prime = @(angle) cos(angle) .* (1 + fn(angle)) ...
                    + sin(angle) .* fn_prime(angle);
                
    % function
    f_objective = @(angle) sin(angle) * (1+ fn(angle)) - r;
    
    %We want a solution that satisfy theta > pi/2!
    if exist('guess', 'var'); theta = guess; else; theta = pi/2 + pi/4; end
    tol_theta = 1e-10;
    n = 1;
    while abs(f_objective(theta)) >= tol_theta && n < 100
       theta = mod(theta - f_objective(theta)/fn_prime(theta), pi);
       n = n + 1;
    end
    %assert(theta > pi/2);

end