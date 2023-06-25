function [unkn, vel] = solve_ODE_unkown(deformation, pressures, dt, previous_conditions, PROBLEM_CONSTANTS)
    if isstruct(deformation); deformation = deformation.deformation_amplitudes; end; if size(deformation, 2) > 1; deformation = deformation'; end
    if isstruct(pressures);   pressures   = pressures.pressure_amplitudes; end; if size(pressures, 2) > 1; pressures = pressures'; end
    nb_harmonics = previous_conditions{end}.nb_harmonics;
    if any(isnan(deformation))|| length(deformation) < nb_harmonics; deformation = zeros(nb_harmonics, 1); calc_vel = true; end
    if any(isnan(pressures))  || length(pressures)   < nb_harmonics;   pressures = zeros(nb_harmonics, 1); calc_vel = false;end
    
    n = length(previous_conditions); % Determines the order of the method
    if n > 2 || n < 1; throw("Hey!"); end

    extract_symbol = @(jj, field) previous_conditions{jj}.(field);
    previous_velocities = zeros(nb_harmonics, n);
    previous_deformation= zeros(nb_harmonics, n);
    for ii = 1:n
        previous_velocities(:, ii) = extract_symbol(ii, 'deformation_velocities');
        previous_deformation(:,ii) = extract_symbol(ii, 'deformation_amplitudes');
    end
    
    if n == 1
        coefs = [-1.0, 1.0];
    elseif n == 2
        rk = dt/previous_conditions{end}.dt;
        ak = (1+2*rk)/(1+rk);
        bk = -(1+rk);
        ck = rk^2/(1+rk);
        coefs = [ck, bk, ak]; 
    end
    om = PROBLEM_CONSTANTS.omegas_frequencies;
    result = (coefs(end)^2/dt + om.^2 * dt) .* deformation + dt * ((1:nb_harmonics)' .* pressures) ...
         + sum(coefs(1:(end-1)) .* (coefs(end) *  previous_deformation/dt + previous_velocities), 2);
     
    if calc_vel
        unkn = -result./(coefs(end)^2/dt + om.^2 * dt);
        vel = (coefs(end) * unkn + sum(coefs(1:(end-1)) .* previous_deformation, 2))/dt;
        vel(1)  = 0; 
        unkn(1) = 0; 
    else
        unkn = -result./(dt * (1:nb_harmonics)');
        vel = nan;
    end
    if size(unkn, 1) > 1; unkn = unkn'; end
    if size(vel,  1) > 1; vel  = vel';  end
end