% Author: Elvis Aguero
% email: elvisavfc65@gmail.com
% Date: January 13th, 2023.

function [recorded_conditions, recorded_times, PROBLEM_CONSTANTS] = solve_motion_v2(varargin)
    
    %
    %Tries to solve the full kinematic match between a dropplet
    % and a solid substrate in vacuum conditions.
    % Arguments:
    % varargin must contain at most three elements, and all of them must be structs whose fields may be left
    % empty or assume nan values to adopt default values, when possible.
    %  - varargin{1} = Physical parameters for the simulation. All dimensions are in cgs units. Can contain the following fields:
    %                                (dim, default) <-- Default shown as ? if variable is required
    %      1) undisturbed_radius     (cm, 1)          = radius of the undisturbed sphere 
    %      2) com_height             (cm, contact)    = Initial height of the COM of the drop. Default is iminent contact
    %      3) initial_velocity       (cm/s, ?)        = Initial velocity of the COM (negative = towards impact)
    %      4) amplitudes             (cm, all zeroes) = Spectral legendre amplitudes for the surface of the drop. 1-indexed
    %      5) amplitudes_vel         (cm, all zeroes) = Spectral legendre velocities for the surface of the drop. 1-indexed
    %      6) pressure_amplitudes    (. , all zeroes) = Initial pressure amplitudes in spectral coordinates. 0-indexed
    %      7) initial_contact_points (adim, 0)        = Initial number of contact points
    %      8) rhoS                   (kg/cm^3, 0.988) = The density of the fluid inside the droplet. Default = water
    %      9) sigmaS                 (., 72.29)       = Surface tension of the fluid inside the droplet. Default = water
    %      10) g                     (cm/s^2, 9.81e+2)= Gravitational constant.
    %      11) nu              (St = cm^2/s, .978e-2) = Kinematic Viscocity
    %  - varargin{2} = Numerical parameters for the simulation
    %      1) harmonics_qtt          (adim, ?)        = Number of spectral amplitudes that describe the motion. 
    %      2) angular_sampling (adim, harmonics_qtt+1)= Number of angles that describe the shape of the drop.
    %      3) simulation_time        (s, inf)         = Maximum simulation time allowed. 
    %                                                   (Inf = simulate roughly until contact has ended)
    %      4) version                (int,  3)        = version for the system of equations to be solved.
    %           v1 = Nonlinear exact integration on contact area (Default)
    %           v2 = Nonlinear approximated integration on whole sphere
    %           v3 = Linearised version of v2 (Only first non constant pressure harmonic contributing)
    %      5) order                  (int, 1)         = order of the finite
    %                                                   differences to be solved
    %  - varargin{3} = Other options. 
    %      1) live_plotting          (bool, true)     = whether or not to plot real-time results (more consuming)
    %      2) debug_flag             (bool, false)    = Verbose real-time info for the simulation (experimental feature)
    %      3) folder                 (string, ") .    = Folder whre the
    %                                 results will be stored. Default is current directory
    %      4) prefix                 (sttring, '')    = Prefix to be put in
    %      the name of the output files
    %      5) save_results           (bool, true)     = Whether to export or not the results
    %      6) optimize_for_bounce    (bool, false)    = Whether to exit early if bounce has ended
    %      7) saving_frequency       (s, nan)         = Determines the frequency at which we save the output results
    
    %% Defining default Arguments
    
    default_options = struct('live_plotting', false, 'debug_flag', false, ...
            'folder', fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), '2_output'), ...
            'prefix', '', 'save_results', true, 'optimize_for_bounce', false, ...
            'saving_frequency', nan);
    default_numerical = struct('simulation_time', inf, 'harmonics_qtt', 20, ...
        'angular_sampling', nan, 'version', 3, 'order', 1);
    
    default_physical = struct('undisturbed_radius', .05, 'initial_height', inf, ...
        'initial_velocity', -10, 'initial_amplitudes', nan, ...
        'amplitudes_velocities', nan, 'pressure_amplitudes', nan, ...
        'initial_contact_points', 0, 'rhoS', 0.988, 'sigmaS', 72.20, ...
        'g', 9.81e+2, 'nu', .978e-2);

    %% Handling default arguments. All units are in cgs.
    if nargin >= 3
        % if isstruct(varargin{3}) == false; error('Option values is not a struct'); end

        % Overriding default values
        A = fieldnames(varargin{3});
        for ii = 1:length(A)
            default_options.(A{ii}) = varargin{3}.(A{ii});
        end
    end
    if nargin >= 2
        
        if isstruct(varargin{2}) == false; error('Numerical values is not a struct'); end
        % overriding default values
        A = fieldnames(varargin{2});
        for ii = 1:length(A)
            default_numerical.(A{ii}) = varargin{2}.(A{ii});
        end
    end
    
    if nargin >= 1
        if isstruct(varargin{1}) == false; error('Numerical values is not a struct'); end
        % Overriding default values
        A = fieldnames(varargin{1});
        for ii = 1:length(A)
            default_physical.(A{ii}) = varargin{1}.(A{ii});
        end
    end
    

    if ~isfinite(default_numerical.harmonics_qtt)
        default_numerical.harmonics_qtt = 20; 
    end
    if ~isfinite(default_numerical.angular_sampling)
        default_numerical.angular_sampling = default_numerical.harmonics_qtt + 1; 
    end
    harmonics_qtt = default_numerical.harmonics_qtt; 
    if isnan(default_physical.initial_amplitudes); default_physical.initial_amplitudes = zeros(1, harmonics_qtt); end
    if isnan(default_physical.amplitudes_velocities); default_physical.amplitudes_velocities = zeros(1, harmonics_qtt); end
    if isnan(default_physical.pressure_amplitudes); default_physical.pressure_amplitudes = zeros(1, harmonics_qtt+1); end
    
    %% Setting default values to variables for conciseness
    undisturbed_radius = default_physical.undisturbed_radius; 
    initial_height = default_physical.initial_height; 
    initial_velocity=default_physical.initial_velocity; 
    initial_amplitudes=default_physical.initial_amplitudes; 
    pressure_amplitudes = default_physical.pressure_amplitudes; 
    initial_contact_points = default_physical.initial_contact_points; 
    rhoS = default_physical.rhoS; 
    sigmaS= default_physical.sigmaS; 
    harmonics_qtt = default_numerical.harmonics_qtt; 
    simulation_time = default_numerical.simulation_time; 
    version=default_numerical.version;
    angular_sampling = default_numerical.angular_sampling;
    g = default_physical.g;
    nu = default_physical.nu;
    debug_flag = default_options.debug_flag;
    live_plotting = default_options.live_plotting;
    prefix = default_options.prefix;
    

    % Dimensionless Units
    length_unit = undisturbed_radius;
    time_unit = sqrt(rhoS * length_unit^3 / sigmaS); %undisturbed_radius/velocity_unit; % Temporal dimensionless number
    velocity_unit = length_unit/time_unit; % abs(initial_velocity);    
    pressure_unit = rhoS * velocity_unit^2; 
    froude_nb   = length_unit/(g*time_unit.^2);
    weber_nb    = rhoS * undisturbed_radius.^3/ (sigmaS*time_unit.^2); % Weber's number of the dropplet
    Oh = nu * sqrt(rhoS/(sigmaS*undisturbed_radius));% Ohnesorge number

    mS = rhoS * length_unit^3;
    mass_unit = mS;
    
    % Definition of angles that correspond to samplings
    f = @(n)  sqrt(n .* (n+2) .* (n-1));
    omegas_frequencies = f(1:harmonics_qtt)';
    %f2 = @(n) (1 - n ./ ((angular_sampling+n):-1:(1+n))) * pi * (angular_sampling + n)/angular_sampling;
    syms x;
    
    % we choose the angles to be the zeros of the last legendre Polynomial (SOUTH POLE BASED)
    theta_vector = [pi; acos(double(vpasolve(legendreP(angular_sampling-1, x))))]'; clear x;

    
    %% Initial conditions
    % Set dropplet's sphere height initial conditions
    get_initial_height = @(amplitudes) 1 - sum(arrayfun(@(idx) amplitudes(idx) * (-1.0)^(idx), 1:length(amplitudes)));
    
    if length(initial_amplitudes) ~= harmonics_qtt
        initial_amplitudes = zeros(1, harmonics_qtt);
    else
        initial_amplitudes = initial_amplitudes/length_unit;
    end
    if initial_height == Inf || isnan(initial_height)
        initial_height = get_initial_height(initial_amplitudes);
    else
        assert(initial_height >= get_initial_height(initial_amplitudes));
        initial_height = initial_height/length_unit;
    end
    
    initial_velocity_adim = initial_velocity/velocity_unit;

    %if length(initial_mplitude_velocities) ~= harmonics_qtt
        initial_mplitude_velocities = zeros(1, harmonics_qtt);
    %else
    %    initial_mplitude_velocities = initial_mplitude_velocities/velocity_unit;
    %end

    initial_pressure_coefficients = pressure_amplitudes / pressure_unit; % Just to emphasize the units of these coefficients.
    contact_points = 0;
    
    % Define the time step so that the highest frequency has N steps
    N = 8;
    max_dt = (2*pi/(sqrt(harmonics_qtt*(harmonics_qtt+2)*(harmonics_qtt-1))*N));
    %max_dt = round(time_unit/(N * harmonics_qtt^(3/2)), 1, 'significant')/time_unit; %max_dt = 0.0034359491980026735;
    dt = max_dt; 
    
    initial_time = 0;
    current_time = initial_time/time_unit;
    if simulation_time == inf
        % Just to save some values in the matrix
        final_time = 20e-3/time_unit; %min(50000*dt, 10e-3/time_unit);
    else
        final_time = simulation_time/time_unit;
    end
    if isnan(default_options.saving_frequency)
        saving_frequency = 1e-5/time_unit; % One shot every  0.01 ms
    else
        saving_frequency = default_options.saving_frequency;
    end
    saving_frequency = max(saving_frequency, max_dt);
    counter_time = saving_frequency;
    current_index = 2;%  This integer points to the next available index in variables that are going to 
                      %  export data (index 1 is for initial conditions)
    maximum_index = ceil(final_time/saving_frequency) + 5;%ceil((final_time - initial_time)/dt) + 4;
    %anumber_of_extra_indexes = 0;

    %grow_dt = false;%  THis variable controls how fast dt can grow
    %jjj = 0; myflag = true;%  Indexes to keep track how small is dt compared to max_dt
    flips = 0;
      
    legendre_matrix = precompute_integrals(theta_vector, harmonics_qtt);
    function_to_minimize = eval(sprintf('@function_to_minimize_v%d', default_numerical.version));
    JacobianCalculator = eval(sprintf('@JacobianCalculator_v%d', default_numerical.version));
    % Constants of the problem formulation
    PROBLEM_CONSTANTS = struct("froude_nb", froude_nb, "weber_nb", weber_nb, ...
        "Oh", Oh, ...
        "version", version, ...
        "nb_harmonics", harmonics_qtt, ...
        "omegas_frequencies", omegas_frequencies, ...
        "angles_qtt", angular_sampling, ... % number of angles 
        "pressure_unit", pressure_unit, ...
        'precomputed_integrals', legendre_matrix, ...
        "theta_vector", theta_vector, ... % set of fixed angle vectors
        "function_to_minimize", function_to_minimize, ... % v1 = fully nonlinear integration on disk, v2 = nonlinear with spherical approximation, v3 = linearised version of v2
        "jacobian_calculator", JacobianCalculator, ... % has to be the same version as function to minimize
        "DEBUG_FLAG", debug_flag); %true = plot and save video
    mat_inverses = struct();

    
    current_conditions = ProblemConditions_v2( ...
        harmonics_qtt, ...
        initial_amplitudes, ...
        initial_mplitude_velocities, ...
        initial_pressure_coefficients, ...
        current_time, ...
        dt(end), ...
        initial_height, ...
        initial_velocity_adim, initial_contact_points); % Last argument is contact radius
 
    previous_conditions = {current_conditions, current_conditions}; 
    % 1-st order set
    previous_conditions = previous_conditions(end);

   %  Preallocate variables that will be exported (All of them have units!)
   recorded_conditions =cell(maximum_index, 1); 
   give_dimensions_v2 = @(X) ProblemConditions_v2( ...
       X.nb_harmonics, ...
        X.deformation_amplitudes * length_unit, ...
        X.deformation_velocities * velocity_unit, ...
        X.pressure_amplitudes * (mass_unit * length_unit / (time_unit^2 * length_unit^2)), ...
        X.current_time * time_unit, ...
        X.dt * time_unit, ...
        X.center_of_mass * length_unit, ...
        X.center_of_mass_velocity * velocity_unit, ...
        X.contact_points); 

     recorded_conditions{1} = give_dimensions_v2(previous_conditions{end});
     recorded_times = zeros(maximum_index, 1); recorded_times(1) =  current_time * time_unit;

    indexes_to_save = zeros(maximum_index, 1); indexes_to_save(1) = 1;
    current_to_save = 2;
    
    file_path = fullfile(default_options.folder, sprintf('Version v%d (rhoS=%.2e, sigmaS=%.2e, R=%.2e)', version, rhoS, sigmaS, undisturbed_radius)); %fullfile(sprintf("../2_output/%s/", default_options.folder));
    if exist(file_path, 'dir') ~= 7 % CHeck if folder simulations exists
        mkdir(file_path); % If not, create it
    end    
    
    datestring = datestr(datetime(), 30);
    if PROBLEM_CONSTANTS.DEBUG_FLAG == true
        file_name2 = fullfile(file_path, sprintf('%s.mp4', datestring));
        vidObj = VideoWriter(file_name2,'MPEG-4');
        set(vidObj,'FrameRate',10)
        open(vidObj);
        
        progress_bar = 0.1;
    end

    
    % Define advancing step function: 
    %   v4 = classic newton method
    %   v5 = Newton method with best value estimator
    advance_one_step = @get_next_step_v5;
    %% Starting main loop
    init = clock;
    
    try
        while ( current_time < final_time) 
            % First, we try to solve with the same number of contact points
            probableNextConditions = cell(5, 1);
            errortan = Inf * ones(1, 5);
            recalculate = false;
            
            [probableNextConditions{3}, errortan(3), mat_inverses] = advance_one_step(previous_conditions, dt(end), ...
                contact_points, PROBLEM_CONSTANTS, mat_inverses);
            [probableNextConditions{4}, errortan(4), mat_inverses] = advance_one_step(previous_conditions, dt(end), ...
                contact_points+1, PROBLEM_CONSTANTS, mat_inverses);
            [probableNextConditions{2}, errortan(2), mat_inverses] = advance_one_step(previous_conditions, dt(end), ...
                contact_points-1, PROBLEM_CONSTANTS, mat_inverses);
                
            if (abs(errortan(3)) > abs(errortan(4)) || abs(errortan(3)) > abs(errortan(2)))
                if abs(errortan(4)) <= abs(errortan(2))
                    %Now lets check with one more point to be sure
                    [~, errortan(5), mat_inverses] = advance_one_step(previous_conditions, dt(end), ...
            contact_points+2, PROBLEM_CONSTANTS, mat_inverses);
    
                    if abs(errortan(4)) < abs(errortan(5))
                        %Accept new data
                        
                        current_conditions = probableNextConditions{4};
                        contact_points = contact_points + 1;
                    else
                        %time step too big
                        recalculate = true;
                    end
                else
                    %now lets check if errortan is good enough with one point
                    %less
                    [~, errortan(1), mat_inverses] = advance_one_step(previous_conditions, dt(end), ...
            contact_points-2, PROBLEM_CONSTANTS, mat_inverses);
    
    
                    if abs(errortan(2)) < abs(errortan(1))
                        %Accept new data
                        current_conditions = probableNextConditions{2};
                        contact_points = contact_points - 1;
                    else
                        recalculate = true;
                    end
                end %End of (errortan(4) < errortan(2))
    
            else %the same number of contact points is best    
                if errortan(3) == Inf % ALl errors are infinity
                    recalculate = true;
                else
                    %Accept new data
                    current_conditions = probableNextConditions{3};
                    
                end
            end %
            
            if recalculate == true
                dt = [dt(1:(end-1)), dt(end)/2, dt(end)/2];
                fprintf("Se dividio dt por 2, dt=%e U = %g, modes = %g, version = %d \n", dt(end), initial_velocity, harmonics_qtt, version);
                % Refine time step in index notation 
                %iii = iii + 1; jjj = 2 * jjj;
                
                if dt(end) * time_unit < 1e-12 % Time step this small is not physically meaningful
                    fprintf("Time step too small (%e). U = %g, modes = %g, version = %d \n", dt(end), initial_velocity, harmonics_qtt, version);
                    error("Time step too small (%e). U = %g, modes = %g, version = %d \n", dt(end), initial_velocity, harmonics_qtt, version);
                end
                % If one hour has elapsed, close this simulation
                if etime(clock, init) > 60 * 120
                    error("Too much time has elapsed U = %g, modes = %g, version = %d \n", initial_velocity, harmonics_qtt, version);
                end
            else
                % Progressively increase order of method until desired
                % length of previous conditions have been attained
                if default_numerical.order > length(previous_conditions) && current_time/dt(end) >= 10
                    previous_conditions = {previous_conditions{1:end} current_conditions};
                else
                    previous_conditions = {previous_conditions{2:end} current_conditions};
                end

                current_time = current_time + dt(end); %jjj = jjj + 1;
                dt = dt(1:max(1, end-1)); if size(dt, 1) == 1; dt = max_dt; end
                % if mod(jjj, 2) == 0 && grow_dt == true
                %     jjj = floor(jjj/2); 
                %     iii = iii - 1;
                %     % Increase time step
                %     dt = 2 * dt;
                %     % Decrease the number of time you can make dt bigger
                %     grow_dt = false;
                % end
                
                % Only store data if you have passed the counter time
                % (multiples of saving_frequency variable above)
                if current_time >= counter_time
                    % TODO: % Stored data
                    recorded_conditions{current_index} = give_dimensions_v2(current_conditions);
                    recorded_times(current_index) = current_time * time_unit;
                    current_index = current_index + 1; % Point to the next space in memory 
                    
                    indexes_to_save(current_to_save) = current_index - 1;
                    current_to_save = current_to_save + 1;
                    % % If we are in a multiple of max_dt, reset indexes
                    % if dt(end) == max_dt
                    % 
                    %     indexes_to_save(current_to_save) = current_index - 1;
                    %     current_to_save = current_to_save + 1;
                    % else
                    %     number_of_extra_indexes = number_of_extra_indexes + 1;
                    % end
                    counter_time = counter_time + saving_frequency;

                    if live_plotting == true 
                        % Do some live plotting here
        
                        plot_title = sprintf(" t = %-8.5f (ms), Contact points = %-2g deg, \n v_k = %-8.5f cm/s, z_k = %-8.5f cm\n", ...
                           1e+3 * current_time * time_unit, current_conditions.contact_points, ...
                                current_conditions.center_of_mass_velocity * velocity_unit, ...
                                current_conditions.center_of_mass* length_unit);
                        h = plot_condition(1, current_conditions, 1.25, plot_title, PROBLEM_CONSTANTS);
                        
                        if PROBLEM_CONSTANTS.DEBUG_FLAG == true
                            currFrame = getframe(h);
                            writeVideo(vidObj,currFrame);
                        end
        
                    else
                        if PROBLEM_CONSTANTS.DEBUG_FLAG == true
                            if current_time/final_time >= progress_bar
                                fprintf("Done %.0f%% of the simulation.\n", current_time/final_time*100);
                                progress_bar = progress_bar + 0.1;
                            end
                        end
                    end   
                end

                if current_index > 2 && ...
                        recorded_conditions{current_index-1}.center_of_mass_velocity * recorded_conditions{current_index-2}.center_of_mass_velocity <= 0
                    flips = flips + 1;
                end
                velocity_down = (current_conditions.center_of_mass_velocity < 0 && current_time > 10*max_dt);
                if (simulation_time == inf &&  contact_points==0 && ...
                        (velocity_down || current_conditions.center_of_mass > 1.1)) || ...
                        (default_options.optimize_for_bounce == true && velocity_down == true && flips > 0)
                    final_time = current_time*1.1;
                    simulation_time = 1e+6; % So as not to enter to this if ever again
                    flips = -inf;
                    if PROBLEM_CONSTANTS.DEBUG_FLAG==true
                        fprintf("Changed final time. Current progress: %.0f%%\n", ...
                        current_time/final_time*100);
                    end
                end
                if simulation_time == inf && current_time >= 0.999*final_time
                    warning("Simulation time is infinite but simulation apparently wont finish")
                    fprintf("Too much time has elapsed U = %g, modes = %g, version = %d \n", initial_velocity, harmonics_qtt, version);
                    break;
                end
    
                
            end
    
        end
        file_name = fullfile(file_path, sprintf("%s%s.mat", ...
           prefix, datestring));
    catch me
        file_name = fullfile(file_path, sprintf("%s%s (Errored).mat", ...
            prefix, datestring));
        PROBLEM_CONSTANTS.error = me;
    end
    
    if PROBLEM_CONSTANTS.DEBUG_FLAG == true
        close(vidObj);
    end
    
    
    % Post processing
    indexes_to_save = indexes_to_save(1:(current_to_save-1));
    recorded_conditions = recorded_conditions(indexes_to_save); recorded_times = recorded_times(indexes_to_save);
    for ii = 1:10
        if exist(file_name, 'file')
            file_name = replace(file_name, ".mat", sprintf("-%d.mat", randi(9)));
        else
            break;
        end
    end
    if default_options.save_results == true
    save(file_name, 'recorded_conditions', 'recorded_times', 'length_unit', ...
        'velocity_unit', 'pressure_unit', 'time_unit', 'froude_nb', 'weber_nb', 'PROBLEM_CONSTANTS', ...
        'default_numerical', 'default_physical'); 
    end
    
end
