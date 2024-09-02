function [unkn, vel] = solve_ODE_unkown(deformation, pressures, dt, previous_conditions, PROBLEM_CONSTANTS)

    % This function will solve the discrete ordinary differential equation
    % on the modes of pressures and amplitudes, whichever is missing
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


    %omegas = PROBLEM_CONSTANTS.omegas_frequencies;
    harmonics = (1:nb_harmonics)';
    amplitudes_coefficients = harmonics .* (harmonics+2) .* (harmonics-1); % amplitude term
    pressure_coeffs = harmonics; % Pressure term
    Oh = PROBLEM_CONSTANTS.Oh;
    vel_coeffs = PROBLEM_CONSTANTS.Oh * 2 * (2*harmonics + 1) .* (harmonics-1); % Viscocity term
    result = (coefs(end)/dt*(coefs(end) + dt * Oh*vel_coeffs) + amplitudes_coefficients * dt) .* deformation + ...
        dt * (pressure_coeffs .* pressures) ...
         + sum(coefs(1:(end-1)) .* ((coefs(end) + dt * Oh*vel_coeffs) .*  previous_deformation/dt + previous_velocities), 2);
     
    if calc_vel
        unkn = -result./(coefs(end)/dt*(coefs(end) + dt * Oh*vel_coeffs) + amplitudes_coefficients * dt);
        vel = (coefs(end) * unkn + sum(coefs(1:(end-1)) .* previous_deformation, 2))/dt;
        vel(1)  = 0; 
        unkn(1) = 0;
    else
        unkn = -result./(dt * harmonics);
        vel = nan;
    end
    if size(unkn, 1) > 1; unkn = unkn'; end
    if size(vel,  1) > 1; vel  = vel';  end
end