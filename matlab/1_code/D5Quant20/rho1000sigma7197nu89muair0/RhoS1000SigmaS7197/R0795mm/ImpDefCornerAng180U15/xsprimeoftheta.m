function xsp = xsprimeoftheta(theta,A2,A3)

c = cos(theta);
s = sin(theta);
c2 = cos(2*theta);
s2 = sin(2*theta);
c4 = cos(4*theta);

xsp = (1+A2/4)*c-A3*c2/4-3*A2*s.*s2/2+3*A2*c.*c2/4+5*A3*c4/4;