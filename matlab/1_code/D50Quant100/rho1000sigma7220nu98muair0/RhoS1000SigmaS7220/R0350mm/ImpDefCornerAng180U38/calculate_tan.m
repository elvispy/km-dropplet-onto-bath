function tann = calculate_tan(distance_to_axis, amplitudes)


    if isstruct(amplitudes); amplitudes = amplitudes.deformation_amplitudes; end
    angle = theta_from_cylindrical(distance_to_axis, amplitudes);
    zeta = zeta_generator(amplitudes);
    if size(amplitudes, 2) > 1; amplitudes = amplitudes'; end
    der = @(theta) sum(amplitudes .* collectdnPl(length(amplitudes), cos(theta)), 1);

    dzdr = @(theta)  (-sin(theta) .* (1 + zeta(theta)) - cos(theta) .* sin(theta) .* der(theta)) ./ ...
        (cos(theta) .* (1 + zeta(theta)) - sin(theta).^2 .* der(theta));
        
    tann =  tan(dzdr(angle));



    % Calculates the tangent of the sphere at distance distance_to_axis(ii)
%     fn = zeta_genera
%     
%     M = length(distance_to_axis);
%     thetas = theta_from_cylindrical(distance_to_axis, amplitudes);
% %     for ii =1:(M-2)
% %         if ii > 1
% %             thetas(ii) = theta_from_cylindrical(distance_to_axis(ii), fn, fn_prime, thetas(ii-1));
% %         else
% %             thetas(ii) = theta_from_cylindrical(distance_to_axis(ii), fn, fn_prime, pi);
% %         end
% %     end
%     RR = arrayfun(fn, thetas);
%     RRprime = arrayfun(fn_prime, thetas); 
%     
%     tann = double((-sin(thetas) .* (1 + RR) + cos(thetas).* RRprime) ./ ...
%             (cos(thetas) .* (1 + RR) + sin(thetas) .* RRprime));
end