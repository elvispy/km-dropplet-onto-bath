function theta = thetaofxs(xs,A2,A3,theta)

tol = 1E-10;

c = cos(theta);
s = sin(theta);
c2 = cos(2*theta);
s2 = sin(2*theta);
s4 = sin(4*theta);
c4 = cos(4*theta);
f = (1+A2/4)*s-A3*s2/8+3*A2*s.*c2/4+5*A3*s4/16-xs;
fprime = (1+A2/4)*c-A3*c2/4-3*A2*s.*s2/2+3*A2*c.*c2/4+5*A3*c4/4;


while abs(f)>tol
    theta = mod(abs(theta-f/fprime),pi);
    c = cos(theta);
    s = sin(theta);
    c2 = cos(2*theta);
    s2 = sin(2*theta);
    s4 = sin(4*theta);
    c4 = cos(4*theta);
    f = (1+A2/4)*s-A3*s2/8+3*A2*s.*c2/4+5*A3*s4/16-xs;
    fprime = (1+A2/4)*c-A3*c2/4-3*A2*s.*s2/2+3*A2*c.*c2/4+5*A3*c4/4;
end
%syms x;
%assert(abs(theta_from_cylindrical(xs, [0, A2, A3], diff(legendreP(1:3, x)))-theta) < 1e-10 || abs(theta_from_cylindrical(xs, [0, A2, A3], diff(legendreP(1:3, x)))+theta-pi) < 1e-10, sprintf("r=%g, A2=%g, A3=%g", xs, A2, A3));
end

function theta = theta_from_cylindrical(r, amplitudes, legendrePprime)
    % Calculates the theta that corresponds to a radius by Newton-Raphson
    % Where fn(theta) = sin(theta) * (R + zeta)
    N = length(amplitudes); syms x;
    
    % Derivative of the function i'm trying to find the zero of
    fn_prime = @(angle) cos(angle) .* (1 + sum(amplitudes .* legendreP(1:N, cos(angle)))) ...
                    - sin(angle).^2 .* double(sum(amplitudes .* subs(legendrePprime, x, cos(angle))));
                
    % function
    fn = @(angle) sin(angle) .* (1 + sum(legendreP(1:N, cos(angle)) .* amplitudes)) - r;
    
    %We want a solution that satisfy theta > pi/2!
    theta = pi/2 + pi/4;
    tol_theta = 1e-10;
    n = 1;
    while abs(fn(theta)) >= tol_theta && n < 100
       theta = mod(theta - fn(theta)/fn_prime(theta), pi);
       n = n + 1;
    end
    %assert(theta > pi/2);

end
