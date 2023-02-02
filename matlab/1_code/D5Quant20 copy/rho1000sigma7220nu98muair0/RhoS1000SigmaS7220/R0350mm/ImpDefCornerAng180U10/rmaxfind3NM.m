function maxRad = rmaxfind3NM(A2,A3)

theta = thetaMax(-pi,A2,A3);

maxRad = sin(theta)+1.5*A2*cos(theta)^2*sin(theta)-.5*A2*sin(theta)...
    +2.5*A3*cos(theta)^3*sin(theta)-1.5*A3*cos(theta)*sin(theta);

