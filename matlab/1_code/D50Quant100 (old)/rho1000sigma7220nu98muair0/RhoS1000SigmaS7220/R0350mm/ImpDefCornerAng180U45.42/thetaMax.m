function theta = thetaMax(A2,A3)

tol = 1E-10;

theta = .9*pi;
c = cos(theta);
s = sin(theta);
c2 = cos(2*theta);
s2 = sin(2*theta);
s4 = sin(4*theta);
c4 = cos(4*theta);

xsprime = (1+A2/4)*c-A3*c2/4-3*A2*s.*s2/2+3*A2*c.*c2/4+5*A3*c4/4; % First derivative
xssecond = -(1+A2/4)*s+A3*s2/2-3*A2*c.*s2-15*A2*s.*c2/4-5*A3*s4; % Second derivative
nn = 1;
while abs(xsprime)>tol && nn < 100% Newton's method
    theta = mod(abs(theta-xsprime/xssecond),pi);
    c = cos(theta);
    s = sin(theta);
    c2 = cos(2*theta);
    s2 = sin(2*theta);
    c4 = cos(4*theta);
    xsprime = (1+A2/4)*c-A3*c2/4-3*A2*s.*s2/2+3*A2*c.*c2/4+5*A3*c4/4;
    xssecond = -(1+A2/4)*s+A3*s2/2-3*A2*c.*s2-15*A2*s.*c2/4-5*A3*s4;
    nn = nn + 1;
end

%fprintf("theta_old = %.15f, thetanew = %.15f", theta, theta_max_radius([0, A2, A3]));
end

function theta = theta_max_radius(amplitudes)
    % Calculates the theta that corresponds to the maximum radius possible
    % Using newton raphson with the equation fn'(theta) = 0
    % where fn(theta) = sin(theta) (R + zeta)
    
    N = length(amplitudes); syms ang x;
    
    % function
    fn =  sin(ang)*(1 +  sum(amplitudes .* legendreP(1:N, cos(ang))));
  
    % Derivative of the funtion
    fn_prime = diff(fn, ang); fn_prime_anonymous = @(angle) double(subs(fn_prime, ang, angle));
    
    % second derivative of the function
    fn_primeprime = diff(fn_prime, ang); fn_primeprime_anonymous = @(angle) double(subs(fn_primeprime, ang, angle));
    
    %We want a solution that satisfy theta > pi/2!
    theta = pi/2 + pi/4;
    tol_theta = 1e-10;
    n = 1;
    while abs(fn_prime_anonymous(theta)) >= tol_theta && n < 100
       theta = mod(theta - fn_prime_anonymous(theta)/fn_primeprime_anonymous(theta), pi);
       n = n + 1;
    end
    %disp(n);
    assert(theta > pi/2-pi/6);

end