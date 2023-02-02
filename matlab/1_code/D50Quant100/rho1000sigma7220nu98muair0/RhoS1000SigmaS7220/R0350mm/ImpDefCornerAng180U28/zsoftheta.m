function zs = zsoftheta(theta,A2,A3)

c = cos(theta);
c2 = cos(2*theta);
c4 = cos(4*theta);

zs = 3*A3/16+(1+A2/4)*c+A3*c2/2+3*A2*c2.*c/4+5*A3*c4/16;
end
function zs = zs_from_spherical(theta, amplitudes) % check pointwise function
    zs = 1;
    s = cos(theta);
    for ii = 1:length(amplitudes)
        zs = zs + amplitudes(ii) * legendreP(ii, s);
    end
    zs = zs .* s;
end