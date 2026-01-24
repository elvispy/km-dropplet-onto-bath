% SOLVE_MOTION  Simulate drop/bath impact using spherical-harmonic dynamics.
%
%   solve_motion(U0, ~, N, tolP, wd, debug_flag, tend_override)
%
%   This routine advances the coupled motion of a deformable droplet
%   impacting a bath modeling surface with N spherical-harmonic modes.
%
%   Refactored 2026-01-24 to use adaptive contact point search (k-logic).

function solve_motion(U0, ~, N, tolP, wd, debug_flag, tend_override)

    lastwarn('', '');
    tstart = tic;

    if nargin < 7 || isempty(tend_override); tend_override = 20; end
    if nargin < 6 || isempty(debug_flag); debug_flag = false; end
    if nargin < 5 || isempty(wd); wd = pwd; end
    
    validateattributes(U0, {'numeric'}, {'scalar', 'real', 'finite'});
    validateattributes(N, {'numeric'}, {'scalar', 'integer', 'positive'});
    validateattributes(tolP, {'numeric'}, {'scalar', 'real', 'positive'});

    Ang = 180; 
    currfold = pwd;
    cleanupObj = onCleanup(@() cd(currfold));

    try
        params = load_simulation_setup(wd, U0, Ang, N, tolP);
    catch ME
        fprintf('Error in load_simulation_setup: %s\n', ME.message);
        error("Working directory not ready to perform simulation");
    end

    % Unpack parameters
    Ro = params.Ro; rhoS = params.rhoS; sigmaS = params.sigmaS;
    rho = params.rho; sigma = params.sigma; nu = params.nu; muair = params.muair;
    g_accel = params.g; D = params.D; quant = params.quant;
    nr = params.nr; dr = params.dr; Delta = params.Delta;
    IntMat = params.IntMat; DTN = params.DTN;
    Ma = params.Ma; Ra = params.Ra;

    % Characteristic Units
    L_unit = Ro; 
    M_unit = rhoS * L_unit^3;
    T = sqrt(rhoS * Ro^3/sigmaS); 
    T_unit = T;
    V_unit = L_unit/T_unit;

    % Dimensionless numbers
    Dr = rhoS/rho; 
    Re = L_unit^2/(nu*T_unit);
    Fr = L_unit/(g_accel * T_unit^2);
    We = rho * L_unit.^3 / (sigma * T_unit^2); 
    WeS  = rhoS*Ro^3/(sigmaS * T_unit^2); 
    Westar = rhoS * U0.^2 * Ro / sigmaS; 
    Oh = nu*sqrt(rhoS/(sigmaS*Ro));
    Cang = (Ang/180)*pi; 

    tend = tend_override; 

    % Initial conditions
    t = 0;
    etao = zeros(nr,1); 
    phio = zeros(nr,1); 

    nsteps = ceil(20/(2*pi) * N^(3/2)); 
    dtb = 1/nsteps; 
    
    % Storage
    tvec = t:(dtb):(tend+1); tvecOri = tvec;
    nTime = numel(tvec);
    
    etaOri = zeros(1, nTime);
    z = zeros(1, nTime);
    vz = zeros(1, nTime);
    numl = zeros(1, nTime);
    
    oscillation_amplitudes = zeros(N, nTime); 
    pressure_amplitudes    = zeros(N, nTime); 
    Rv = -ones(1, nTime);
    oscillation_velocities_storage = zeros(N, nTime); 

    nlmax = zeros(1, nTime);
    etas      = zeros(length(etao), nTime); etas(:, 1) = etao;
    phis      = zeros(length(etao), nTime); phis(:, 1) = phio;
    psMatPer = cell(1, nTime);
    psMatPer{1} = zeros(quant+1,1);
    
    indexes_to_save = zeros(nTime, 1);
    indexes_to_save(1) = 1;
    current_to_save = 2;
    savingvarNames = {'z', 'etaOri', 'etas', 'psMatPer', 'vz', 'tvec', ...
        'nlmax', 'numl', 'oscillation_amplitudes', 'pressure_amplitudes', 'Rv'};

    % Frequencies
    f = @(n) sqrt(n.*(n+2).*(n-1)./WeS);
    omegas_frequencies = f(1:N)';

    % Initial State Setup
    amplitudes_old = oscillation_amplitudes(:, 1);
    amplitudes_velocities_old = oscillation_velocities_storage(:, 1);
    B_l_ps_old = zeros(1, N);

    z(1) = -1* zs_from_spherical(pi, oscillation_amplitudes(:, 1));
    vz(1) = -abs(U0/ V_unit); 

    current_conditions = struct("deformation_amplitudes", amplitudes_old, ...
        "deformation_velocities", amplitudes_velocities_old, ...
        "pressure_amplitudes", B_l_ps_old, "dt", dtb, "nb_harmonics", N,  ...
        "current_time", 0, ...
        "center_of_mass", z(1), "center_of_mass_velocity", vz(1), ...
        "nb_contact_points", 0);

    previous_conditions = {current_conditions, current_conditions}; 
    
    % Back-propagate for ODE starter
    dt = dtb;
    previous_conditions{1}.current_time = previous_conditions{2}.current_time - dt;
    previous_conditions{1}.center_of_mass_velocity = ...
        previous_conditions{2}.center_of_mass_velocity + dt/Fr;
    previous_conditions{1}.center_of_mass = ...
        previous_conditions{2}.center_of_mass - previous_conditions{2}.center_of_mass_velocity * dt;

    mode_solution = @(t, idx) current_conditions.deformation_amplitudes(idx) * cos(f(idx) * t) ...
        + current_conditions.deformation_velocities(idx)/(f(idx)+1e-30) * sin(f(idx) * t); 

    for idx = 1:N
        previous_conditions{1}.deformation_amplitudes(idx) = mode_solution(-dt, idx);
        previous_conditions{1}.deformation_velocities(idx) = (mode_solution(0, idx) - mode_solution(-2*dt/1000, idx))/(2*dt/1000);
    end

    PROBLEM_CONSTANTS = struct("froude_nb", Fr, "weber_nb", We, ...
        "nb_harmonics", N, "density_ratio", Dr, "omegas_frequencies", omegas_frequencies, ...
        "spatial_tol", dr, "initial_dt", dtb, "DEBUG_FLAG", debug_flag, "linear_on_theta", true, ...
        "Ra", Ra, "interpolation_number", 10, 'Oh', Oh, 'Westar', Westar);
    
    PROBLEM_CONSTANTS.IntMat = IntMat; PROBLEM_CONSTANTS.Re = Re; 
    PROBLEM_CONSTANTS.Delta = Delta; PROBLEM_CONSTANTS.DTN = DTN;
    PROBLEM_CONSTANTS.Ma = Ma; PROBLEM_CONSTANTS.Cang = Cang;
    PROBLEM_CONSTANTS.nr = nr; PROBLEM_CONSTANTS.dr = dr;
    PROBLEM_CONSTANTS.Fr = Fr; PROBLEM_CONSTANTS.We = We;

    fprintf("Starting simulation on %s\n", pwd);

    tentative_index = 0;
    going_back = 0;
    
    % Tracking variables locally
    curr_eta = etas(:,1);
    curr_phi = phis(:,1);
    curr_z = z(1);
    curr_vz = vz(1);
    curr_numl = 0;
    curr_osc_amps = oscillation_amplitudes(:,1);
    curr_ps = psMatPer{1};
    
    % State struct to pass to attempt_step (simplifies arguments)
    state = struct('eta', curr_eta, 'phi', curr_phi, 'z', curr_z, 'vz', curr_vz, ...
                   'numl', curr_numl, 'osc_amps', curr_osc_amps, ...
                   'ps', curr_ps, 'prev_cond', {previous_conditions});

    %% Main Loop
    try
        while t < tend
            tentative_index = tentative_index + 1;
            ensure_storage_capacity(tentative_index + 1);

            t = tvec(tentative_index+1);
            dt = t - tvec(tentative_index);
            
            % if PROBLEM_CONSTANTS.DEBUG_FLAG
               % fprintf("Time: %0.4g, dt: %0.3e, tentative_idx: %d\n", t, dt, tentative_index); 
            % end
            
            recalculate = true;
            loop_limit = 0;
            
            while recalculate
                recalculate = false;
                loop_limit = loop_limit + 1;
                if loop_limit > 50
                     error("Too many retries required for a single step.");
                end

                 [~, warnId] = lastwarn();
                 if strcmp(warnId, 'MATLAB:nearlySingularMatrix') || dt*T_unit < 1e-10
                      recalculate = true; 
                      lastwarn('', '');
                      going_back = going_back + 1;
                      if PROBLEM_CONSTANTS.DEBUG_FLAG
                          fprintf("Warning detected (small dt or singular). Going back count: %d\n", going_back);
                      end
                      if going_back > 50
                          error("Went back too many times (dt too small or singular). Stopping execution.");
                      end
                 else
                     % New wrapper returns accepted result or sets recalculate=true
                     best_res = attempt_step_wrapper(state.numl, state.numl, dt, state, PROBLEM_CONSTANTS, N, tolP);
                     recalculate = best_res.recalculate;
                 end
                 
                 if recalculate
                     % Reduce dt
                     % We must be careful not to keep modifying tvec infinitely for the same index without effect
                     % insert_midpoint reduces dt for the current step
                     tvec = insert_midpoint(tvec, tentative_index);
                     tentative_index = tentative_index - 1; 
                     break; 
                 end
                 
                 if ~isempty(best_res)
                     state = update_state(state, best_res, dt, t);
                     
                     idx_save = tentative_index + 1;
                     numl(idx_save) = state.numl;
                     z(idx_save) = state.z;
                     vz(idx_save) = state.vz;
                     etaOri(idx_save) = state.eta(1);
                     etas(:, idx_save) = state.eta;
                     phis(:, idx_save) = state.phi;
                     psMatPer{idx_save} = state.ps;
                     nlmax(idx_save) = best_res.nlmaxTent;
                     Rv(idx_save) = best_res.RvTent;
                     oscillation_amplitudes(:, idx_save) = state.osc_amps;
                     pressure_amplitudes(:, idx_save) = state.prev_cond{2}.pressure_amplitudes;
                     
                     if t >= tvecOri(current_to_save) || t < 10*dtb
                         indexes_to_save(current_to_save) = idx_save;
                         current_to_save = current_to_save + 1;
                     end
                     
                     if (state.z > 1.5 && state.numl == 0) || (state.vz < 0 && vz(tentative_index) >= 0 && t > 50 * dtb)
                        tend = t;
                     end
                     
                     if PROBLEM_CONSTANTS.DEBUG_FLAG
                         debug_plot(state, best_res, t, PROBLEM_CONSTANTS);
                     end
                 end
            end % while recalculate
            
         end % while t < tend
    catch ME
        save_result_safe("errored_");
        fprintf("Simulation Errored: %s\n", ME.message);
        rethrow(ME);
    end

    save_result_safe("");
    
    simul_time = toc(tstart);
    save('ProblemConditions.mat', "T", "N", "U0", "Ang", "Re", "Fr", "We", ...
         "WeS", "Cang", "tend", "nsteps", "dtb", "L_unit", "T_unit", "M_unit", ...
         "PROBLEM_CONSTANTS", "simul_time");
    fprintf("Finished simulation. Time: %0.2f min\n", simul_time/60);

    % --- Nested Functions ---
    
    function res_out = attempt_step_wrapper(prev_k, k, dt, st, PC, N, tolP)
        % Generalized state machine for contact point search using k-logic
        
        res_out = struct('err', inf, 'recalculate', true); 
        
        % Helper to find nlprev based on history to capture hysteresis/advancing state
        function np = find_nlprev(target_k)
            % Search backwards from tentative_index for a value different from target_k
            if tentative_index == 0
                np = 0; return;
            end
            
            hist_vec = numl(tentative_index:-1:1);
            co = find(hist_vec ~= target_k, 1);
            
            if isempty(co)
                np = target_k;
            else
                np = numl(tentative_index - co + 1);
            end
        end

        % Helper to run one trial
        function r = try_k(candidate_k)
             if candidate_k < 0
                 r = struct('err', inf); return;
             end
             prev_k_hist = find_nlprev(candidate_k);
             r = attempt_step(prev_k_hist, candidate_k, dt, st, PC, N, tolP); 
             if PC.DEBUG_FLAG
                 fid = fopen('debug_solve_motion.log', 'a');
                 fprintf(fid, "Time: %g, dt: %g, try_k(%d) [nlprev=%d] -> err: %g\n", st.prev_cond{2}.current_time, dt, candidate_k, prev_k_hist, r.err);
                 fclose(fid);
             end
        end

        % 1. Evaluate Current State (k) and Neighbors (k-1, k+1)
        r0 = try_k(k);
        r_minus = try_k(k-1);
        r_plus = try_k(k+1);

        err0 = abs(r0.err);
        errM = abs(r_minus.err);
        errP = abs(r_plus.err);

        % Optimization: If current state is flight (k=0) and good enough, stay.
        if k == 0 && err0 < 0.5
            res_out = r0; res_out.recalculate = false; return;
        end

        % 2. Compare and Decide
        if err0 <= errP && err0 <= errM
             % Current k is best among neighbors.
             if err0 == inf
                 res_out.recalculate = true; return;
             else
                 res_out = r0; res_out.recalculate = false; return;
             end
             
        elseif errP < err0 && errP <= errM
             % k+1 is better. Check stability (k+2).
             r_plus2 = try_k(k+2);
             if abs(r_plus2.err) > errP
                  % k+1 is a local minimum (better than k and k+2)
                  res_out = r_plus; res_out.recalculate = false; return;
             else
                  % k+2 is even better? Step size too large or rapid change. Recalculate.
                  res_out.recalculate = true; return;
             end

        else % errM < err0 && errM < errP
             % k-1 is better. Check stability (k-2).
             r_minus2 = try_k(k-2);
             if abs(r_minus2.err) > errM
                  % k-1 is a local minimum
                  res_out = r_minus; res_out.recalculate = false; return;
             else
                  % k-2 is even better. Recalculate.
                  res_out.recalculate = true; return;
             end
        end
    end

    function res = attempt_step(prev_k, k, dt, st, PC, N, tolP)
        res = struct();
        
        amps_curr = st.osc_amps;
        RmaxOld = r_from_spherical(maximum_contact_radius(amps_curr), amps_curr);
        nlmax_pred = max(floor(RmaxOld/PC.dr)+1, st.numl);
        thetaVec_pred = theta_from_cylindrical(PC.dr*(0:(nlmax_pred-1)), amps_curr);
        
        B_l_ps_tent = project_pressure_amplitudes(st.ps, st.numl, thetaVec_pred, PC.dr, amps_curr, N, PC);
        [amplitudes_tent, ~] = solve_ODE_unkown(nan, B_l_ps_tent, dt, st.prev_cond, PC);
        
        [nlmaxTent, thetaVec, zs, RvTent, angleDropMP] = compute_contact_geometry(amplitudes_tent, amps_curr, amplitudes_tent, PC.dr, PC.nr);
        
        if k > nlmaxTent
             res.err = inf;
             return;
        end
        
        psRes = [];
        
        curr_eta = st.eta; curr_phi = st.phi; curr_z = st.z; curr_vz = st.vz;
        
        % Iterate Pressure
        errorP = 1; iter = 0;
        
        while errorP >= tolP && iter < 100
             iter = iter + 1;
             
             if k == 0
                 % solveDD0 expects: solveDD0(dt,z,vz,etao,phio,nr,Re,Delta,DTN,Fr,We,zs,Rv)
                 % PC struct has fields, but we passed them individually incorrectly before if PC.Fr or PC.We were missing
                 [etaTent, phiTent, zTent, vzTent, errT] = solveDD0(dt, curr_z, curr_vz, curr_eta, curr_phi, PC.nr, PC.Re, PC.Delta, PC.DTN, PC.Fr, PC.We, zs, RvTent);
                 psRes = zeros(nlmaxTent, 1);
                 % Enforce threshold for non-contact
                 if abs(errT) > 0.5; errT = inf; end
             else
                  [etaTent, phiTent, zTent, vzTent, psRes, errT] = solvenDDCusp(prev_k, k, dt, curr_z, curr_vz, curr_eta, curr_phi, PC.nr, PC.dr, PC.Re, PC.Delta, PC.DTN, PC.Fr, PC.We, 0, zs, PC.IntMat(k,:), angleDropMP, PC.Cang, PC.density_ratio, RvTent);
             end
             
             B_l_ps_new = project_pressure_amplitudes(psRes, k, thetaVec, PC.dr, amps_curr, N, PC);
             [amplitudes_new, velocities_new] = solve_ODE_unkown(nan, B_l_ps_new, dt, st.prev_cond, PC);
             
             errorP = norm(amplitudes_tent - amplitudes_new)/norm(amplitudes_tent);
             amplitudes_tent = amplitudes_new;
        end
       
        res.err = abs(errT);
        res.k = k; % Store the number of contact points used
        res.eta = etaTent; res.phi = phiTent; res.z = zTent; res.vz = vzTent;
        res.ps = psRes; res.amps = amplitudes_new; res.vels = velocities_new;
        res.B_l_ps = B_l_ps_new; res.nlmaxTent = nlmaxTent; res.RvTent = RvTent; res.zs = zs;
    end
    
    function new_st = update_state(st, res, dt, t)
        new_st = st;
        new_st.eta = res.eta; new_st.phi = res.phi;
        new_st.z = res.z; new_st.vz = res.vz;
        new_st.ps = res.ps; new_st.osc_amps = res.amps;
        
        % Trust the decision logic's choice of k
        new_st.numl = res.k;
        
        new_st.prev_cond{1} = new_st.prev_cond{2};
        new_st.prev_cond{2} = struct("deformation_amplitudes", res.amps, ...
            "deformation_velocities", res.vels, ...
            "dt", dt, "nb_harmonics", N, "pressure_amplitudes", res.B_l_ps, ...
            "current_time", new_st.prev_cond{1}.current_time + dt, ...
            "center_of_mass", res.z, "center_of_mass_velocity", res.vz, ...
            "nb_contact_points", new_st.numl);
    end
    
     function ensure_storage_capacity(requiredIndex)
        currentLen = size(z, 2);
        if requiredIndex <= currentLen; return; end
        newLen = max(requiredIndex, currentLen + ceil(0.2 * currentLen) + 16);
        etaOri(1, newLen) = 0; z(1, newLen) = 0; vz(1, newLen) = 0; numl(1, newLen) = 0; nlmax(1, newLen) = 0;
        indexes_to_save(newLen, 1) = 0;
        oscillation_amplitudes(:, newLen) = 0; pressure_amplitudes(:, newLen) = 0;
        etas(:, newLen) = 0; phis(:, newLen) = 0;
        if numel(Rv) < newLen; Rv(1, newLen) = -1; Rv((currentLen+1):newLen) = -1; end
        if numel(psMatPer) < newLen; psMatPer{newLen} = []; end
        oscillation_velocities_storage(:, newLen) = 0; 
    end

    function save_result_safe(prefix)
         indexes = indexes_to_save(1:(current_to_save-1));
         s = struct(); s.z = z(indexes); s.etaOri = etaOri(indexes); s.etas = etas(:, indexes);
         s.psMatPer = psMatPer(indexes); s.vz = vz(indexes); s.tvec = tvec(indexes);
         s.nlmax = nlmax(indexes); s.numl = numl(indexes);
         s.oscillation_amplitudes = oscillation_amplitudes(:, indexes);
         s.pressure_amplitudes = pressure_amplitudes(:, indexes);
         s.Rv = Rv(indexes);
         for n = fieldnames(s)'
             % save -struct expects a scalar structure containing the fields to be saved.
             % We want to save a file named e.g. 'z.mat' containing a variable named 'z'.
             % So we constructing s_sub with field 'z'.
             % This logic seems correct IF save is called correctly.
             % The error was: Argument to -STRUCT must be the name of a scalar structure variable.
             % This happens if 's_sub' is not seen as a scalar structure? It is.
             % Or maybe because we construct it inline? No.
             % WAIT: save(..., '-struct', 's_sub') requires 's_sub' to be a VARIABLE NAME string in the workspace?
             % Let's use the explicit save(filename, varname) pattern, it's safer.
             % But we need to rename the variable in the workspace to match n{1}.
             % eval is dirty but works for dynamic naming in save.
             % Or struct path.
             
             % Cleaner:
             temp_struct = struct();
             temp_struct.(n{1}) = s.(n{1});
             save(sprintf('%s%s.mat', prefix, n{1}), '-struct', 'temp_struct');
         end
    end
    
    function debug_plot(st, res, t, PC)
         if ~isempty(get(groot, 'CurrentFigure'))
             set(gcf, 'Visible', 'off');
         else
             figure('Visible', 'off');
         end
         
         xs = PC.dr*(0:res.nlmaxTent-1);
         zsplot = res.zs(1:res.nlmaxTent)+res.RvTent+st.z;
         plot([-fliplr(xs(2:end)),xs],[flipud(zsplot(2:end));zsplot],'k','Linewidth',2);
         hold on
         width = min(PC.nr, 200); xplot = PC.dr * (0:PC.nr-1);
         plot([-fliplr(xplot(2:width)),xplot(1:width)],[flipud(st.eta(2:width));st.eta(1:width)],'LineWidth',2);
         hold off
         axis equal; title(sprintf('t = %0.3f', t)); grid on; drawnow;
    end
end

function params = load_simulation_setup(wd, U0, Ang, N, tolP)
    if isstring(wd); wd = char(wd); end
    validateattributes(wd, {'char'}, {'nonempty'});
    run_dir_hint = wd; impact_dir = fileparts(run_dir_hint); r0_dir = fileparts(impact_dir);
    solid_dir = fileparts(r0_dir); fluid_dir = fileparts(solid_dir); domain_dir = fileparts(fluid_dir);
    Ro = load_matvar(fullfile(r0_dir, 'Ro.mat'), 'Ro');
    rhoS = load_matvar(fullfile(solid_dir, 'rhoS.mat'), 'rhoS'); sigmaS = load_matvar(fullfile(solid_dir, 'sigmaS.mat'), 'sigmaS');
    rho = load_matvar(fullfile(fluid_dir, 'rho.mat'), 'rho'); sigma = load_matvar(fullfile(fluid_dir, 'sigma.mat'), 'sigma');
    nu = load_matvar(fullfile(fluid_dir, 'nu.mat'), 'nu'); muair = load_matvar(fullfile(fluid_dir, 'muair.mat'), 'muair');
    g_accel = load_matvar(fullfile(fluid_dir, 'g.mat'), 'g');
    D = load_matvar(fullfile(domain_dir, 'D.mat'), 'D'); quant = load_matvar(fullfile(domain_dir, 'quant.mat'), 'quant');
    nr = load_matvar(fullfile(domain_dir, 'nr.mat'), 'nr'); dr = load_matvar(fullfile(domain_dir, 'dr.mat'), 'dr');
    Delta = load_matvar(fullfile(domain_dir, 'Delta.mat'), 'Delta'); IntMat = load_matvar(fullfile(domain_dir, 'IntMat.mat'), 'IntMat');
    dtn_path = fullfile(domain_dir, sprintf('DTNnew345nr%dD%drefp10.mat', nr, D)); DTN = load_matvar(dtn_path, 'DTNnew345');
    xplot = dr * (0:nr-1); save(fullfile(domain_dir, 'xplot.mat'), 'xplot');

    fluid_name = ['rho', num2str(1000*rho), 'sigma', num2str(round(100*sigma)), 'nu', num2str(round(10000*nu)), 'muair', num2str(muair)];
    solid_name = ['RhoS', num2str(rhoS*1000), 'SigmaS', num2str(round(100*sigmaS))];
    r0_name = sprintf('R0%gmm', Ro*10000); impact_name = sprintf("ImpDefCornerAng%gU%.4g", Ang, U0);
    run_name = sprintf('N=%dtol=%0.2e', N, tolP);
    run_dir = fullfile(domain_dir, fluid_name, solid_name, r0_name, impact_name, run_name);
    if exist(run_dir, 'dir') ~= 7; error('Run directory not found: %s', run_dir); end
    Ma = load_matvar(fullfile(domain_dir, fluid_name, solid_name, 'Ma.mat'), 'Ma');
    Ra = load_matvar(fullfile(domain_dir, fluid_name, solid_name, 'Ra.mat'), 'Ra');
    cd(run_dir);
    params = struct("Ro", Ro, "rhoS", rhoS, "sigmaS", sigmaS, "rho", rho, "sigma", sigma, "nu", nu, ...
        "muair", muair, "g", g_accel, "D", D, "quant", quant, "nr", nr, "dr", dr, ...
        "Delta", Delta, "IntMat", IntMat, "DTN", DTN, "Ma", Ma, "Ra", Ra, "run_dir", run_dir);
end

function value = load_matvar(path, var_name)
    if exist(path, 'file') ~= 2; error('Missing file: %s', path); end
    loaded = load(path, var_name);
    if ~isfield(loaded, var_name); error('Missing variable %s in %s', var_name, path); end
    value = loaded.(var_name);
end

function B_l_ps = project_pressure_amplitudes(ps, nb_contact_points, thetaVec, dr, amplitudes_for_r, N, PROBLEM_CONSTANTS)
    if isempty(ps) || nb_contact_points < 1 || norm(ps, 1) == 0; B_l_ps = zeros(1, N); return; end
    nb_contact_points = min(nb_contact_points, numel(ps));
    if PROBLEM_CONSTANTS.linear_on_theta == true
        if nb_contact_points > 1 && numel(thetaVec) >= nb_contact_points
            contactAngle = (1.5 * thetaVec(nb_contact_points) - 0.5 * thetaVec(nb_contact_points-1));
        elseif numel(thetaVec) >= 2
            contactAngle = (thetaVec(2) + thetaVec(1)) / 2;
        else
            contactAngle = thetaVec(1);
        end
        angles_qtt = max(2, (nb_contact_points + 1) * PROBLEM_CONSTANTS.interpolation_number);
        angles = linspace(contactAngle, pi, angles_qtt);
        radii = r_from_spherical(angles, amplitudes_for_r);
        f = @(r) interp1(dr*(0:nb_contact_points), [ps(1:nb_contact_points)', 0], r, 'linear', 0);
        values = f(radii);
        B_l_ps = custom_project_amplitudes(angles, values, N, NaN, NaN);
    else
        % Fallback
        if nb_contact_points >= numel(thetaVec)
             if numel(thetaVec) < 2; angles = [thetaVec(:)', pi]; else; angles = [thetaVec(1:nb_contact_points), (2 * thetaVec(nb_contact_points) - thetaVec(nb_contact_points-1))]; end
        else; angles = thetaVec(1:(nb_contact_points+1)); end
        f = @(thetas) interp1(angles, [ps(1:nb_contact_points)', 0], thetas, 'linear', 0);
        endpoints = [angles(end), angles(1)];
        B_l_ps = project_amplitudes(f, N, endpoints, PROBLEM_CONSTANTS, true);
    end
end

function [nlmaxTent, thetaVec, zs, RvTent, angleDropMP] = compute_contact_geometry(amplitudes_for_radius, amplitudes_for_theta, amplitudes_for_shape, dr, nr)
    RmaxTent = r_from_spherical(maximum_contact_radius(amplitudes_for_radius), amplitudes_for_radius);
    nlmaxTent = max(1, floor(RmaxTent/dr) + 1);
    thetaVec = theta_from_cylindrical(dr * (0:(nlmaxTent-1)), amplitudes_for_theta);
    RvTent = zs_from_spherical(pi, amplitudes_for_shape);
    zs = zeros(nr, 1);
    zs(1:nlmaxTent) = zs_from_spherical(thetaVec, amplitudes_for_shape)' - RvTent;
    zs((nlmaxTent+1):nr) = Inf;
    tanDrop = calculate_tan(dr * (1:nlmaxTent) - dr/2, amplitudes_for_shape)';
    angleDropMP = atan(tanDrop(1:nlmaxTent));
end

function tvec = insert_midpoint(tvec, idx)
    tvec = [tvec(1:idx), (tvec(idx) + tvec(idx+1))/2, tvec(idx+1:end)];
end
