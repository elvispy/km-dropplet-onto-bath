function zs = rs_from_spherical(theta, zeta)
    if size(theta, 1) > 1; theta = theta'; end
    zs = sin(theta) .* (1 + zeta(theta));
    
%     s = cos(theta);
%     for ii = 1:length(amplitudes)
%         zs = zs + amplitudes(ii) .* legendrep(ii, s);
%     end
%     zs = zs .* sin(theta);
end