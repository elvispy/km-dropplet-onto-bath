function ysp = ysprimeoftheta(theta,A2,A3)

c = cos(theta);
s = sin(theta);
c2 = cos(2*theta);
s2 = sin(2*theta);
s4 = sin(4*theta);

ysp = -(1+A2/4)*s-A3*s2-3*A2*s.*c2/4-3*A2*c.*s2/2-5*A3*s4/4;