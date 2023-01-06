function xs = xsoftheta(theta,A2,A3)

ss = sin(theta);
c2 = cos(2*theta);
s2 = sin(2*theta);
s4 = sin(4*theta);

xs = (1+A2/4)*ss-A3*s2/8+3*A2*c2.*ss/4+5*A3*s4/16;

end
function zs = rs_from_spherical(theta, amplitudes)
    zs = 1;
    s = cos(theta);
    for ii = 1:length(amplitudes)
        zs = zs + amplitudes(ii) .* legendreP(ii, s);
    end
    zs = zs .* sin(theta);
end
