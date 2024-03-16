function thetaMax = xs3(theta,A2,A3)

tol = 1E-6;

xsprime = (1-3.5*A2)*cos(theta)-1.5*A3*cos(2*theta)-15*A3*sin(2*theta)^2/8 ...
    +4.5*A2*cos(theta)^3+2.5*A3*cos(theta)^4;
xssecond = -11*A2*sin(theta)+sin(2*theta)-15*A3*sin(4*theta)/4+...
    27*a2*sin(theta)^3/2-5*A3*cos(theta)^2*sin(2*theta);

while abs(xsprime)>tol
    theta = theta-xsprime/xssecond;
    xsprime = (1-3.5*A2)*cos(theta)-1.5*A3*cos(2*theta)-15*A3*sin(2*theta)^2/8 ...
    +4.5*A2*cos(theta)^3+2.5*A3*cos(theta)^4;
xssecond = -11*A2*sin(theta)+sin(2*theta)-15*A3*sin(4*theta)/4+...
    27*a2*sin(theta)^3/2-5*A3*cos(theta)^2*sin(2*theta);
end
thetaMax = theta;