function zs = zsoftheta(theta,A2,A3)

c = cos(theta);
c2 = cos(2*theta);
c4 = cos(4*theta);

% ys = 3*A3/16+(1+.25*A2)*c+c2.*(3*A2*c/4+A3/2)-5*A3*c4/16;
zs = 3*A3/16+(1+A2/4)*c+A3*c2/2+3*A2*c2.*c/4+5*A3*c4/16;