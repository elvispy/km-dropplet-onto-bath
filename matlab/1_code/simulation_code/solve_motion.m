% SOLVE_MOTION  Simulate drop?bath impact using spherical-harmonic dynamics.
%
%   solve_motion(U0, ~, N, tolP, wd, debug_flag)
%
%   This routine advances the coupled motion of a deformable droplet
%   impacting a bath, modeling the free surface with N spherical-harmonic
%   modes and iterating a pressure/contact problem at each time step until
%   convergence. It adaptively refines the time grid when tangency/pressure
%   checks fail, saves time histories of key fields, and writes a summary of
%   run conditions on exit.
%
%   INPUTS
%   ------
%   U0          Impact speed of the droplet [cm/s in CGS]. Converted to
%               dimensionless units internally using Ro and sigmaS, rhoS.
%
%   ~           Unused placeholder (kept for interface consistency).
%
%   N           Integer. Number of spherical-harmonic modes retained.
%
%   tolP        Positive scalar. Convergence tolerance for the
%               inner pressure?shape fixed-point loop (relative change of
%               deformation amplitudes).
%
%   wd          Char/string. Path to the working directory where the run?s
%               parameter subfolders and cache files live. The function
%               `cd`s into `wd` and navigates relative to it to locate inputs.
%
%   debug_flag  Logical. If true, prints per-step diagnostics and live plots
%               of the surface and contact region.
%
%   EXPECTED INPUT FILES (MAT) & DIRECTORY LAYOUT
%   ---------------------------------------------
%   The routine expects a parameterized folder tree rooted at `wd`'s parent
%   directories and loads the following files (all in CGS unless noted):
%
%     Ro.mat          : Ro  (sphere/drop radius)
%     rhoS.mat        : rhoS (sphere density)
%     sigmaS.mat      : sigmaS (sphere surface tension)
%     rho.mat         : rho  (fluid density)
%     sigma.mat       : sigma (fluid surface tension)
%     nu.mat          : nu   (kinematic viscosity)
%     muair.mat       : muair (air viscosity; used for path selection)
%     g.mat           : g    (gravity)
%
%     D.mat           : D       (domain diameter in units of Ro)
%     quant.mat       : quant   (# of dr?s per Ro for output sampling)
%     nr.mat          : nr      (# radial grid points)
%     dr.mat          : dr      (radial grid spacing)
%     Delta.mat       : Delta   (time-integration parameter)
%     IntMat.mat      : IntMat  (contact-integral stencils)
%     DTNnew345nr%D%drefp10.mat : DTNnew345 (Dirichlet-to-Neumann operator)
%
%     Under folder names keyed by fluid/solid parameters:
%       Ma.mat        : Ma   (dimensionless mass)
%       Ra.mat        : Ra   (density ratio for coupling)
%
%   The code constructs subfolder names such as:
%     ./rho{1000*rho}sigma{round(100*sigma)}nu{round(10000*nu)}muair{muair}
%     /RhoS{rhoS*1000}SigmaS{round(100*sigmaS)}/R0{Ro*10000 in mm}
%     /ImpDefCornerAng{Ang}U{U0}/N={N}tol={tolP}
%   and will error if these are missing.
%
%   MODEL & NUMERICS (high level)
%   -----------------------------
%   ? Non-dimensionalization uses Ro (length), T = sqrt(rhoS*Ro^3/sigmaS)
%     (time), and corresponding V_unit = Ro/T.
%   ? Governing deformation coordinates are the first N SH modes with
%     natural frequencies ?_l = sqrt(l(l+2)(l-1)/WeS).
%   ? Main outer loop steps time with a baseline dtb ? O(N^{-3/2}) and
%     adaptively inserts midpoints upon:
%       - tangency errors at the free surface,
%       - nearly singular linear systems (lastwarn:
%         'MATLAB:nearlySingularMatrix'),
%       - nonconvergence of the contact/pressure iteration.
%   ? At each time step:
%       1) Predict deformation via ODE step with tentative pressure
%          amplitudes B_l.
%       2) Build/update contact geometry (number of pressed nodes `numl`,
%          max radius, local slope/tangent).
%       3) Solve one of several contact cases (0?4 points, interior/boundary)
%          via `solveDD0` / `solvenDDCusp` to get ?, ?, z, v_z and a new
%          discrete pressure `ps`.
%       4) Project `ps` to SH pressure amplitudes B_l and re-integrate the
%          deformation ODE.
%       5) Iterate until ?A_new ? A_old?/?A_old? < tolP or refine dt.
%
%   OUTPUTS / SIDE EFFECTS
%   ----------------------
%   This function does not return variables; it **saves** per-run artifacts
%   into the current parameter folder:
%
%     z.mat                     : center of mass height vs. time
%     etaOri.mat                : surface height below south pole
%     etas.mat                  : full surface elevation snapshots (nr × T)
%     psMatPer.mat              : accepted pressure distributions (cell)
%     vz.mat                    : center of mass velocity vs. time
%     tvec.mat                  : (possibly refined) time vector
%     nlmax.mat                 : max potential contact nodes per step
%     numl.mat                  : accepted contact nodes per step
%     oscillation_amplitudes.mat: SH deformation amplitudes (N × T)
%     pressure_amplitudes.mat   : SH pressure amplitudes (N × T)
%     Rv.mat                    : south-pole height of the undeformed cap
%
%   Additionally:
%     ProblemConditions.mat : summary of nondimensional numbers, units,
%                             N, U0, Ang, dtb, runtime, and struct
%                             PROBLEM_CONSTANTS used in the simulation.
%     error_logU0=... .mat  : saved if an exception occurs (with ME).
%     Partial results with prefix 'errored_' are saved on errors.
%
%   TERMINATION
%   -----------
%   The run ends when either:
%     ? t reaches an internal cap (default earliest t_end ? 9 non-diml),
%     ? the drop is in free flight with z > 1.5 and no contact, or
%     ? v_z changes sign to negative after sufficient elapsed time.
%
%   WARNINGS & ERROR HANDLING
%   -------------------------
%   ? Nearly singular systems or dt*T_unit < 1e-10 trigger step rollback,
%     time-grid pruning, and midpoint insertion. After 50 rollbacks, the
%     run aborts with an error.
%   ? If the required folder structure/files are missing, the function
%     throws: "Working directory not ready to perform simulation".
%
%   UNITS
%   -----
%   Inputs are in CGS where applicable; outputs are saved in dimensionless
%   form unless obviously dimensional (e.g., raw parameters mirrored from
%   the .mat files).
%
%   DEPENDENCIES (called internally)
%   --------------------------------
%     solveDD0, solvenDDCusp, solve_ODE_unkown,
%     r_from_spherical, zs_from_spherical, theta_from_cylindrical,
%     calculate_tan, custom_project_amplitudes / project_amplitudes.
%
%   EXAMPLE
%   -------
%     % Prepare parameter folders/files (see above), then run:
%     U0 = 38;     % cm/s
%     N  = 20;     % modes
%     tolP = 1e-4; % pressure/deformation tolerance
%     wd = pwd;    % inside the run folder chain
%     debug = true;
%     solve_motion(U0, [], N, tolP, wd, debug);
%
%   NOTES
%   -----
%   ? The function mutates the working directory (cd) but restores the
%     original folder on exit.
%   ? Ang (contact angle) is currently fixed to 180° inside the code.
%   ? If a preexisting 'z.mat' is found, the code is set up to warn about
%     overwriting (message currently commented).


function solve_motion(U0, ~, N, tolP, wd, debug_flag)

%Reset warning
lastwarn('', '');

tstart = tic;
%data in cgs

if exist('z.mat', 'file') == 2
   % error("Exporting data is going to be overwritten. Please re-allocate files to avoid loss of data");
end


%U0 = 38; %impact velocity in cm/s (unit of velocity for the problem)
Ang = 180; %contact angle to be imposed

currfold = pwd;
try
    cd(wd);
    cd ..
    cd ..
    load('Ro.mat','Ro')%Sphere's radius in CGS

    cd ..
    load('rhoS.mat','rhoS')%Sphere density
    load('sigmaS.mat', 'sigmaS') %Sphere's surface tension

    cd ..
    load('rho.mat','rho')
    load('sigma.mat','sigma')
    load('nu.mat','nu')
    load('muair.mat', 'muair')
    load('g.mat','g') %gravitational constant

    cd ..
    load('D.mat', 'D')%Domain diameter in units of droplet radii
    load('quant.mat' ,'quant')%number of dr's contained in an undeformed dropelt radius
    load('nr.mat','nr')
    load('dr.mat','dr')
    load('Delta.mat','Delta')
    load('IntMat.mat','IntMat')
    load(sprintf('DTNnew345nr%dD%drefp10.mat', nr, D),'DTNnew345')
    DTN = DTNnew345;
    clear DTNnew345
    xplot = dr*(0:nr-1); save('xplot.mat','xplot')%I might remove or relocate this

    cd(['rho',num2str(1000*rho),'sigma',num2str(round(100*sigma)),'nu',num2str(round(10000*nu)),'muair',num2str(muair)])

    cd(['RhoS',num2str(rhoS*1000),'SigmaS',num2str(round(100*sigmaS))])
    load('Ma.mat','Ma')%Dimensionless mass of sphere
    load('Ra.mat','Ra')%Density ratio

    cd(sprintf('R0%gmm',Ro*10000))

    cd(sprintf("ImpDefCornerAng%gU%.4g", Ang, U0))
    cd(sprintf('N=%dtol=%0.2e', N, tolP));
    
catch
    error("Working directory not ready to perform simulation");
end

% #--- 
%N = 20; % Number of harmonics contributing to the oscillation
% #---0

%Characteristic Units
L_unit = Ro; 
M_unit = rhoS * L_unit^3;
T = sqrt(rhoS * Ro^3/sigmaS); % Characteristic time
T_unit = T;
V_unit = L_unit/T_unit;

%Dimensionless numbers for equations
Dr = rhoS/rho; %Sr = sigmaS/sigma;
Re = L_unit^2/(nu*T_unit);
Fr = L_unit/(g * T_unit^2);
We = rho * L_unit.^3 / (sigma * T_unit^2); 
WeS  = rhoS*Ro^3/(sigmaS * T_unit^2); %This is for the bath/dropplet interaction.
Oh = nu*sqrt(rhoS/(sigmaS*Ro));

Cang = (Ang/180)*pi; %contact angle to be imposed

%Physical parameter
tend = 20; %Earliest possible end of simulation in characteristic units

%Inintial conditions for the fluid
t = 0;
etao = zeros(nr,1); %initial surface elevation
phio = zeros(nr,1); %initial surface potential

%Numerical Simulation parameters
nsteps = ceil(20/(2*pi) * N^(3/2)); %minimum number of timesteps in one unit of time
%nsteps = nsteps + 100 - mod(nsteps, 100);
dtb = 1/nsteps; %basic timestep (gets halved as needed over impacts)
steps = ceil((tend-t)/dtb); %estimated minimum number of timesteps

%Zeroing result storing variables
etaOri = zeros(1,steps+1);%height of the surface below the south pole
z = zeros(1,steps+1);%height of the centre of mass
vz = zeros(1,steps+1);%speed of the centre of mass
numl = zeros(1,steps+1);%number of pressed mesh points at each time step
tvec = t:(dtb):(tend+1); tvecOri = tvec;%vector of times assuming no refinement has happened
%plus some extra time just in case the simulation needs to run longer

dt = tvec(2) - tvec(1); indexes_to_save = zeros(steps + 1, 1);
current_to_save = 2; indexes_to_save(1) = 1;
oscillation_amplitudes = zeros(N, steps + 1); % Variable to store
pressure_amplitudes    = zeros(N, steps + 1); % Pressure amplitudes
Rv = -ones(1, steps+1);
% the time dependent amplitude of all the SH
oscillation_velocities = zeros(N, steps+1);
nlmax = zeros(1,steps+1);%Variable to store the number of nodes spanned by the deformed droplet
%zeroing variable that records each part of the sequence of surface states
%etaMatPer = zeros(length(etao),nsteps); 
etas      = zeros(length(etao), steps); etas(:, 1) = etao;
phis      = zeros(length(etao), steps); phis(:, 1) = phio;
%phiMatPer = zeros(length(phio),nsteps);
psMatPer = cell(1,nsteps); % ??? why cell
%Storing initial surface state
%etaMatPer(:,1) = etao;
%phiMatPer(:,1) = phio;
psMatPer{1} = zeros(quant+1,1);

%zeroing the ceiling functions
zs = zeros(nr,1);


%tolP = 1E-4; %error tolerance for  deformation 

%Drop oscillation frequencies
% #--- 
f = @(n) sqrt(n.*(n+2).*(n-1)./WeS);
omegas_frequencies = f(1:N)';

% #---oscillation_amplitudes = zeros(N, steps + 1);
amplitudes_old = oscillation_amplitudes(:, 1);
amplitudes_velocities_old = oscillation_velocities(:, 1);
B_l_ps_old = zeros(1, N);

% # ---
z(1) = -1* zs_from_spherical(pi, oscillation_amplitudes(:, 1));% -1*zsoftheta(pi,A2(1),A3(1)); %height of the centre of mass (CoM) in dimensionless units,

% zsoftheta(pi,A2(1),A3(1)) gives the height of the south pole with
% respect to the CoM, z(1) is chosen so that the drop is just about to touch down
vz(1) = -abs(U0/ V_unit); %Initial velocity of the CoM in dimesionless units


current_conditions = struct("deformation_amplitudes", amplitudes_old, ...
    "deformation_velocities", amplitudes_velocities_old, ...
    "pressure_amplitudes", B_l_ps_old, "dt", dt, "nb_harmonics", N,  ...
    "current_time", 0, ...
    "center_of_mass", z(1), "center_of_mass_velocity", vz(1), ...
    "nb_contact_points", 0);

previous_conditions = {current_conditions, current_conditions}; 

previous_conditions{1}.current_time = previous_conditions{2}.current_time - dt;
previous_conditions{1}.center_of_mass_velocity = ...
    previous_conditions{2}.center_of_mass_velocity + dt/Fr;
previous_conditions{1}.center_of_mass = ...
    previous_conditions{2}.center_of_mass - previous_conditions{2}.center_of_mass_velocity * dt;

g = @(t, idx) current_conditions.deformation_amplitudes(idx) * cos(f(idx) * t) ...
    + current_conditions.deformation_velocities(idx)/(f(idx)+1e-30) * sin(f(idx) * t); 

for idx = 1:N
    previous_conditions{1}.deformation_amplitudes(idx) = g(-dt, idx);
    previous_conditions{1}.deformation_velocities(idx) = (g(0, idx) - g(-2*dt/1000, idx))/(2*dt/1000);
end

tentative_index = 0; %iteration counter
going_back = 0;

errortan = zeros(5,steps+1);%tangency error recorder
ps_accepted = [];

%If there were some initial pressure acting on the surface and sphere I
%would have to change this bit here to reflect the presure distribution


PROBLEM_CONSTANTS = struct("froude_nb", Fr, "weber_nb", We, ...
    "nb_harmonics", N, ...
    "density_ratio", Dr, ...
    "omegas_frequencies", omegas_frequencies, ...
    "spatial_tol", dr, "initial_dt", dtb, ...
    "DEBUG_FLAG", debug_flag, "linear_on_theta", true, ...
    "Ra", Ra, "interpolation_number", 10, ...
    'Oh', Oh);
                            %"pressure_unit", pressure_unit, ...
                            %"CM", 9, ...
                            %"PG", 2, ...
                            %"KILL_OUTSIDE", true, ...
                            %"wigner3j", {precomputed_wigner(harmonics_qtt)}, ...

                            
fprintf("Starting simulation on %s\n", pwd);


% Names of the variables to be stored
savingvarNames = { ...
    getVarName(z), ...
    getVarName(etaOri), ...
    getVarName(etas), ... %getVarName(phis), ...
    getVarName(psMatPer), ...
    getVarName(vz), ...
    getVarName(tvec), ...
    getVarName(nlmax), ...
    getVarName(numl), ...
    getVarName(oscillation_amplitudes), ...
    getVarName(pressure_amplitudes), ...
    getVarName(Rv)...
};

variableValues = cell(size(savingvarNames));


exit = false;
%% Main Loop
try
    while (t<tend) && exit == false

        tentative_index = tentative_index+1;
        t = tvec(tentative_index+1);
        dt = t - tvec(tentative_index);

        if PROBLEM_CONSTANTS.DEBUG_FLAG == true
           fprintf("Outside %0.4g, %0.3e\n", t-dt, dt); 
        end
        %zeroing the tentative solution variables
        etaprob = zeros(nr,5);
        phiprob = zeros(nr,5);
        vzprob  = zeros(1,5);
        zprob   = zeros(1,5);
        errortan(:,tentative_index+1) = 4*ones(5,1);

        psTent = ps_accepted; %Tentative pressure distribution (we start with the previous pressure)

        %Shape at the start of the time step
        %x that corresponds to the max pressed radius     

        RmaxOld = r_from_spherical(maximum_contact_radius(oscillation_amplitudes(:, tentative_index)), oscillation_amplitudes(:, tentative_index));

        %(i.e. where the tangent plane to the droplet is vertical)
        nlmax(tentative_index) = max(floor(RmaxOld/dr)+1, numl(tentative_index));%max number of contact points

        thetaVec = theta_from_cylindrical(dr*(0:(nlmax(tentative_index)-1)), oscillation_amplitudes(:, tentative_index)); % zeros(1,nlmax(jj));%initialising vector of angles of pressed positions

        % Spherical Harmonics modes
        if norm(psTent,1) == 0

            B_l_ps_tent = zeros(1, N);
        else
            nb_contact_points = numl(tentative_index);% nlmax(tentative_index)-find(flipud(psTent),1)+1; %Number of nodes contact points%
            %needs to be integrated against SH modes

            if PROBLEM_CONSTANTS.linear_on_theta == true
                if nb_contact_points > 1
                    contactAngle = (1.5 * thetaVec(nb_contact_points) - 0.5*thetaVec(nb_contact_points-1));
                else
                    contactAngle = (thetaVec(2) + thetaVec(1))/2;
                end
                angles = linspace(contactAngle, ...% + (thetaVec(nb_contact_points) - thetaVec(nb_contact_points-1))/2
                            pi, (nb_contact_points +1) * PROBLEM_CONSTANTS.interpolation_number);
                values = r_from_spherical(angles, oscillation_amplitudes(:, tentative_index));
                f = @(r) interp1(dr*(0:(nb_contact_points)), ...
                    [psTent(1:nb_contact_points)', 0], r, 'linear',  0);
                values = f(values);
                B_l_ps_tent = custom_project_amplitudes(angles, values, N, NaN, NaN);

            else
                % Linear on theta not assumed
                % Defining where the pressure distribution will be
                % projected
                if nb_contact_points == length(thetaVec)
                    angles = [thetaVec(1:(nb_contact_points)), (2 * thetaVec(nb_contact_points) - thetaVec(nb_contact_points - 1))];
                else
                    angles = thetaVec(1:(nb_contact_points+1));
                end
                f = @(thetas) interp1(angles, [psTent(1:nb_contact_points)', 0], thetas, 'linear',  0); 
                endpoints = [angles(end), angles(1)];
                B_l_ps_tent = project_amplitudes(f, N, endpoints, PROBLEM_CONSTANTS, true); 
            end
        end    

        [amplitudes_tent, ~] = solve_ODE_unkown(nan, B_l_ps_tent, dt, ...
            previous_conditions, PROBLEM_CONSTANTS);


        RmaxTent = r_from_spherical(maximum_contact_radius(oscillation_amplitudes(:, tentative_index)), oscillation_amplitudes(:, tentative_index)); 

        nlmaxTent = floor(RmaxTent/dr)+1;
        thetaVec  = theta_from_cylindrical(dr*(0:(nlmaxTent-1)), oscillation_amplitudes(:, tentative_index)); % zeros(1,nlmaxTent);

        RvTent = zs_from_spherical(pi, amplitudes_tent);
        zs(1:nlmaxTent) = zs_from_spherical(thetaVec, amplitudes_tent)' - RvTent;
        zs((nlmaxTent+1):nr) = Inf;

        tanDrop = calculate_tan( dr * (1:nlmaxTent) - dr/2, amplitudes_tent)';
        angleDropMP(1:(nlmaxTent)) = atan(tanDrop(1:(nlmaxTent)));


        psprob = zeros(nlmaxTent,5);%zeroing the vector of potential pressures
        errorP = 1; %error in the pressure field and amplitude of modes
        reduc = 0; %indicator of whether there was a reduction in the time-step size or not
        ll = 0; % Limiting while loop to 100 iterations
        while abs(errorP)>=tolP && reduc == 0 
            ll = ll + 1;


            if numl(tentative_index) < .5 %i.e. if previously in flight (I need to define this as integer)
                [etaprob(:,3),phiprob(:,3),zprob(3),vzprob(3),errortan(3,tentative_index+1)] = ...
                    solveDD0(dt,z(tentative_index),vz(tentative_index),etas(:, tentative_index), phis(:, tentative_index),nr,Re,Delta,DTN,Fr,We,zs,RvTent);
                if abs(errortan(3,tentative_index+1))<.5
                    numlTent = 0;
                    etaTent = etaprob(:,3);
                    phiTent = phiprob(:,3);
                    psNew = zeros(nlmaxTent,1);
                    zTent = zprob(3);
                    vzTent = vzprob(3);
                else
                    co = find(numl(tentative_index:-1:1)~=1,1);
                    [etaprob(:,4),phiprob(:,4),zprob(4),vzprob(4),psprob(1,4),errortan(4,tentative_index+1)] = ...
                        solvenDDCusp(numl(tentative_index-co+1),1,dt,z(tentative_index),vz(tentative_index),etas(:, tentative_index), phis(:, tentative_index),nr,dr,Re,Delta,DTN,Fr,...
                        We,Ma,zs,IntMat(1,:),angleDropMP,Cang,Dr,RvTent);
                    co = find(numl(tentative_index:-1:1)~=2,1);
                    [~,~,~,~,~,errortan(5,tentative_index+1)] = ...    
                        solvenDDCusp(numl(tentative_index-co+1),2,dt,z(tentative_index),vz(tentative_index),etas(:, tentative_index), phis(:, tentative_index),nr,dr,Re,Delta,DTN,Fr,...
                        We,Ma,zs,IntMat(2,:),angleDropMP,Cang,Dr,RvTent);
                    if abs(errortan(4,tentative_index+1)) < abs(errortan(5,tentative_index+1))
                        numlTent = 1;
                        etaTent = etaprob(:,4);
                        phiTent = phiprob(:,4);
                        psNew = psprob(:,4);
                        zTent = zprob(4);
                        vzTent = vzprob(4);
                    else
                        tvec = [tvec(1:tentative_index),tvec(tentative_index)/2+tvec(tentative_index+1)/2,tvec(tentative_index+1:end)];
                        tentative_index = tentative_index-1;
                        reduc = 1;
                    end
                end
            elseif numl(tentative_index)>.5 && numl(tentative_index)<1.5 % i.e. the last number of contact points was 1
                [etaprob(:,2),phiprob(:,2),zprob(2),vzprob(2),errortan(2,tentative_index+1)] = ...
                    solveDD0(dt,z(tentative_index),vz(tentative_index),etas(:, tentative_index), phis(:, tentative_index),nr,Re,Delta,DTN,Fr,We,zs,RvTent);
                if abs(errortan(2,tentative_index+1))<.5
                    numlTent = 0;
                    etaTent = etaprob(:,2);
                    phiTent = phiprob(:,2);
                    psNew = zeros(nlmaxTent,1);
                    zTent = zprob(2);
                    vzTent = vzprob(2);
                else
                    co = find(numl(tentative_index:-1:1)~=1,1);
                    [etaprob(:,3),phiprob(:,3),zprob(3),vzprob(3),psprob(1,3),errortan(3,tentative_index+1)] = ...                     
                        solvenDDCusp(numl(tentative_index-co+1),1,dt,z(tentative_index),vz(tentative_index),etas(:, tentative_index), phis(:, tentative_index),nr,dr,Re,Delta,DTN,Fr,...
                        We,Ma,zs,IntMat(1,:),angleDropMP,Cang,Dr,RvTent);
                    co = find(numl(tentative_index:-1:1)~=2,1);
                    [etaprob(:,4),phiprob(:,4),zprob(4),vzprob(4),psprob(1:2,4),errortan(4,tentative_index+1)] = ...    
                        solvenDDCusp(numl(tentative_index-co+1),2,dt,z(tentative_index),vz(tentative_index),etas(:, tentative_index), phis(:, tentative_index),nr,dr,Re,Delta,DTN,Fr,...
                        We,Ma,zs,IntMat(2,:),angleDropMP,Cang,Dr,RvTent);
                    if abs(errortan(3,tentative_index+1)) < abs(errortan(4,tentative_index+1))
                        numlTent = 1;
                        etaTent = etaprob(:,3);
                        phiTent = phiprob(:,3);
                        psNew = psprob(:,3);
                        zTent = zprob(3);
                        vzTent = vzprob(3);
                    else
                        co = find(numl(tentative_index:-1:1)~=3,1);
                        [~,~,~,~,~,errortan(5,tentative_index+1)] = ...
                            solvenDDCusp(numl(tentative_index-co+1),3,dt,z(tentative_index),vz(tentative_index),etas(:, tentative_index), phis(:, tentative_index),nr,dr,Re,Delta,DTN,Fr,...
                            We,Ma,zs,IntMat(3,:),angleDropMP,Cang,Dr,RvTent);
                        if abs(errortan(4,tentative_index+1)) < abs(errortan(5,tentative_index+1))
                            numlTent = 2;
                            etaTent = etaprob(:,4);
                            phiTent = phiprob(:,4);
                            psNew = psprob(:,4);
                            zTent = zprob(4);
                            vzTent = vzprob(4);
                        else
                            tvec = [tvec(1:tentative_index),tvec(tentative_index)/2+tvec(tentative_index+1)/2,tvec(tentative_index+1:end)];
                            tentative_index = tentative_index-1; 
                            reduc = 1;
                        end
                    end
                end
            elseif numl(tentative_index) > 1.5 && numl(tentative_index) < 2.5 %i.e. the last contact had two points
                [~,~,~,~,errortan(1,tentative_index+1)] = ...
                    solveDD0(dt,z(tentative_index),vz(tentative_index),etas(:, tentative_index), phis(:, tentative_index),nr,Re,Delta,DTN,Fr,We,zs,RvTent);
                if abs(errortan(1,tentative_index+1))<.5
                    tvec = [tvec(1:tentative_index),tvec(tentative_index)/2+tvec(tentative_index+1)/2,tvec(tentative_index+1:end)];
                    tentative_index = tentative_index-1; 
                    reduc = 1;
                else
                    co = find(numl(tentative_index:-1:1)~=2,1);
                    [etaprob(:,3),phiprob(:,3),zprob(3),vzprob(3),psprob(1:2,3),errortan(3,tentative_index+1)] = ...    
                        solvenDDCusp(numl(tentative_index-co+1),2,dt,z(tentative_index),vz(tentative_index),etas(:, tentative_index), phis(:, tentative_index),nr,dr,Re,Delta,DTN,Fr,...
                        We,Ma,zs,IntMat(2,:),angleDropMP,Cang,Dr,RvTent);
                    co = find(numl(tentative_index:-1:1)~=1,1);
                    [etaprob(:,2),phiprob(:,2),zprob(2),vzprob(2),psprob(1,2),errortan(2,tentative_index+1)] = ...    
                        solvenDDCusp(numl(tentative_index-co+1),1,dt,z(tentative_index),vz(tentative_index),etas(:, tentative_index), phis(:, tentative_index),nr,dr,Re,Delta,DTN,Fr,...
                        We,Ma,zs,IntMat(1,:),angleDropMP,Cang,Dr,RvTent);
                    if abs(errortan(2,tentative_index+1)) < abs(errortan(3,tentative_index+1))
                        numlTent = 1;
                        etaTent = etaprob(:,2);
                        phiTent = phiprob(:,2);
                        psNew = psprob(:,2);
                        zTent = zprob(2);
                        vzTent = vzprob(2);
                    else
                        co = find(numl(tentative_index:-1:1)~=3,1);
                        [etaprob(:,4),phiprob(:,4),zprob(4),vzprob(4),psprob(1:3,4),errortan(4,tentative_index+1)] = ...    
                            solvenDDCusp(numl(tentative_index-co+1),3,dt,z(tentative_index),vz(tentative_index),etas(:, tentative_index), phis(:, tentative_index),nr,dr,Re,Delta,DTN,Fr,...
                            We,Ma,zs,IntMat(3,:),angleDropMP,Cang,Dr,RvTent);
                        if abs(errortan(3,tentative_index+1)) < abs(errortan(4,tentative_index+1))
                            numlTent = 2;
                            etaTent = etaprob(:,3);
                            phiTent = phiprob(:,3);
                            psNew = psprob(:,3);
                            zTent = zprob(3);
                            vzTent = vzprob(3);
                        else
                            co = find(numl(tentative_index:-1:1)~=4,1);
                            [~,~,~,~,~,errortan(5,tentative_index+1)] = ...    
                                solvenDDCusp(numl(tentative_index-co+1),4,dt,z(tentative_index),vz(tentative_index),etas(:, tentative_index), phis(:, tentative_index),nr,dr,Re,Delta,DTN,Fr,...
                                We,Ma,zs,IntMat(4,:),angleDropMP,Cang,Dr,RvTent);
                            if abs(errortan(4,tentative_index+1)) < abs(errortan(5,tentative_index+1))
                                numlTent = 3;
                                etaTent = etaprob(:,4);
                                phiTent = phiprob(:,4);
                                psNew = psprob(:,4);
                                zTent = zprob(4);
                                vzTent = vzprob(4);
                            else
                                tvec = [tvec(1:tentative_index),tvec(tentative_index)/2+tvec(tentative_index+1)/2,tvec(tentative_index+1:end)];
                                tentative_index = tentative_index-1; 
                                reduc = 1;
                            end
                        end
                    end
                end
            elseif numl(tentative_index)>2.5 && numl(tentative_index)<nlmaxTent-1.5 %if the last number of contact points was far from the boundaries
                co = find(numl(tentative_index:-1:1)~=numl(tentative_index),1);
                [etaprob(:,3),phiprob(:,3),zprob(3),vzprob(3),psprob(1:numl(tentative_index),3),errortan(3,tentative_index+1)] = ...    
                    solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index),dt,z(tentative_index),vz(tentative_index),etas(:, tentative_index), phis(:, tentative_index),nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(numl(tentative_index),:),angleDropMP,Cang,Dr,RvTent);
                co = find(numl(tentative_index:-1:1)~=numl(tentative_index)-1,1);
                [etaprob(:,2),phiprob(:,2),zprob(2),vzprob(2),psprob(1:numl(tentative_index)-1,2),errortan(2,tentative_index+1)] = ...    
                    solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index)-1,dt,z(tentative_index),vz(tentative_index),etas(:, tentative_index), phis(:, tentative_index),nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(numl(tentative_index)-1,:),angleDropMP,Cang,Dr,RvTent);
                if abs(errortan(2,tentative_index+1)) < abs(errortan(3,tentative_index+1))
                    co = find(numl(tentative_index:-1:1)~=numl(tentative_index)-2,1);
                    [~,~,~,~,~,errortan(1,tentative_index+1)] = ...
                        solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index)-2,dt,z(tentative_index),vz(tentative_index),etas(:, tentative_index), phis(:, tentative_index),nr,dr,Re,Delta,DTN,Fr,...
                        We,Ma,zs,IntMat(numl(tentative_index)-2,:),angleDropMP,Cang,Dr,RvTent);
                    if abs(errortan(2,tentative_index+1)) < abs(errortan(1,tentative_index+1))
                        numlTent = numl(tentative_index)-1;
                        etaTent = etaprob(:,2);
                        phiTent = phiprob(:,2);
                        psNew = psprob(:,2);
                        zTent = zprob(2);
                        vzTent = vzprob(2);
                    else
                        tvec = [tvec(1:tentative_index),tvec(tentative_index)/2+tvec(tentative_index+1)/2,tvec(tentative_index+1:end)];
                        tentative_index = tentative_index-1; 
                        reduc = 1;
                    end
                else
                    co = find(numl(tentative_index:-1:1)~=numl(tentative_index)+1,1);
                    [etaprob(:,4),phiprob(:,4),zprob(4),vzprob(4),psprob(1:numl(tentative_index)+1,4),errortan(4,tentative_index+1)] = ...    
                        solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index)+1,dt,z(tentative_index),vz(tentative_index),etas(:, tentative_index), phis(:, tentative_index),nr,dr,Re,Delta,DTN,Fr,...
                        We,Ma,zs,IntMat(numl(tentative_index)+1,:),angleDropMP,Cang,Dr,RvTent);
                    if abs(errortan(3,tentative_index+1))<abs(errortan(4,tentative_index+1))
                        numlTent = numl(tentative_index);
                        etaTent = etaprob(:,3);
                        phiTent = phiprob(:,3);
                        psNew = psprob(:,3);
                        zTent = zprob(3);
                        vzTent = vzprob(3);
                    else
                        co = find(numl(tentative_index:-1:1)~=numl(tentative_index)+2,1);%I think I don't need this and I can just replace the first argument of solven by numl(jj)
                        [~,~,~,~,~,errortan(5,tentative_index+1)] = ...
                            solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index)+2,dt,z(tentative_index),vz(tentative_index),etas(:, tentative_index), phis(:, tentative_index),nr,dr,Re,Delta,DTN,Fr,...
                            We,Ma,zs,IntMat(numl(tentative_index)+2,:),angleDropMP,Cang,Dr,RvTent);
                        if abs(errortan(4,tentative_index+1)) < abs(errortan(5,tentative_index+1))
                            numlTent = numl(tentative_index)+1;
                            etaTent = etaprob(:,4);
                            phiTent = phiprob(:,4);
                            psNew = psprob(:,4);
                            zTent = zprob(4);
                            vzTent = vzprob(4);
                        else
                            tvec = [tvec(1:tentative_index),tvec(tentative_index)/2+tvec(tentative_index+1)/2,tvec(tentative_index+1:end)];
                            tentative_index = tentative_index-1; 
                            reduc = 1;
                        end
                    end
                end
            elseif numl(tentative_index) > nlmax(tentative_index)-1.5 && numl(tentative_index) < nlmaxTent-.5 %i.e. if last number of contacted points is nlmax-1
                co = find(numl(tentative_index:-1:1)~=numl(tentative_index),1);
                [etaprob(:,3),phiprob(:,3),zprob(3),vzprob(3),psprob(1:numl(tentative_index),3),errortan(3,tentative_index+1)] = ...    
                    solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index),dt,z(tentative_index),vz(tentative_index),etas(:, tentative_index), phis(:, tentative_index),nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(numl(tentative_index),:),angleDropMP,Cang,Dr,RvTent);
                co = find(numl(tentative_index:-1:1)~=numl(tentative_index)-1,1);
                [etaprob(:,2),phiprob(:,2),zprob(2),vzprob(2),psprob(1:numl(tentative_index)-1,2),errortan(2,tentative_index+1)] = ...    
                    solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index)-1,dt,z(tentative_index),vz(tentative_index),etas(:, tentative_index), phis(:, tentative_index),nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(numl(tentative_index)-1,:),angleDropMP,Cang,Dr,RvTent);
                if abs(errortan(2,tentative_index+1))<abs(errortan(3,tentative_index+1))
                    co = find(numl(tentative_index:-1:1)~=numl(tentative_index)-2,1);
                    [~,~,~,~,~,errortan(1,tentative_index+1)] = ...    
                        solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index)-2,dt,z(tentative_index),vz(tentative_index),etas(:, tentative_index), phis(:, tentative_index),nr,dr,Re,Delta,DTN,Fr,...
                        We,Ma,zs,IntMat(numl(tentative_index)-2,:),angleDropMP,Cang,Dr,RvTent);
                    if abs(errortan(2,tentative_index+1)) < abs(errortan(1,tentative_index+1))
                        numlTent = numl(tentative_index)-1;
                        etaTent = etaprob(:,2);
                        phiTent = phiprob(:,2);
                        psNew = psprob(:,2);
                        zTent = zprob(2);
                        vzTent = vzprob(2);
                    else
                        tvec = [tvec(1:tentative_index),tvec(tentative_index)/2+tvec(tentative_index+1)/2,tvec(tentative_index+1:end)];
                        tentative_index = tentative_index-1; 
                        reduc = 1;
                    end
                else
                    co = find(numl(tentative_index:-1:1)~=numl(tentative_index)+1,1);
                    [etaprob(:,4),phiprob(:,4),zprob(4),vzprob(4),psprob(1:numl(tentative_index)+1,4),errortan(4,tentative_index+1)] = ...    
                        solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index)+1,dt,z(tentative_index),vz(tentative_index),etas(:, tentative_index), phis(:, tentative_index),nr,dr,Re,Delta,DTN,Fr,...
                        We,Ma,zs,IntMat(numl(tentative_index)+1,:),angleDropMP,Cang,Dr,RvTent);
                    if abs(errortan(3,tentative_index+1)) < abs(errortan(4,tentative_index+1))
                        numlTent = numl(tentative_index);
                        etaTent = etaprob(:,3);
                        phiTent = phiprob(:,3);
                        psNew = psprob(:,3);
                        zTent = zprob(3);
                        vzTent = vzprob(3);
                    else
                        numlTent = numl(tentative_index)+1;
                        etaTent = etaprob(:,4);
                        phiTent = phiprob(:,4);
                        psNew = psprob(:,4);
                        zTent = zprob(4);
                        vzTent = vzprob(4);
                    end
                end
            elseif numl(tentative_index) == nlmaxTent %i.e. if last number of contact points was nlmax
                co = find(numl(tentative_index:-1:1)~=numl(tentative_index),1);
                [etaprob(:,3),phiprob(:,3),zprob(3),vzprob(3),psprob(1:numl(tentative_index),3),errortan(3,tentative_index+1)] = ...    
                    solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index),dt,z(tentative_index),vz(tentative_index),etas(:, tentative_index), phis(:, tentative_index),nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(numl(tentative_index),:),angleDropMP,Cang,Dr,RvTent);
                co = find(numl(tentative_index:-1:1)~=numl(tentative_index)-1,1);
                [etaprob(:,2),phiprob(:,2),zprob(2),vzprob(2),psprob(1:numl(tentative_index)-1,2),errortan(2,tentative_index+1)] = ...    
                    solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index)-1,dt,z(tentative_index),vz(tentative_index),etas(:, tentative_index), phis(:, tentative_index),nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(numl(tentative_index)-1,:),angleDropMP,Cang,Dr,RvTent);
                if abs(errortan(2,tentative_index+1)) < abs(errortan(3,tentative_index+1))
                    co = find(numl(tentative_index:-1:1)~=numl(tentative_index)-2,1);
                    [~,~,~,~,~,errortan(1,tentative_index+1)] = ...
                        solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index)-2,dt,z(tentative_index),vz(tentative_index),etas(:, tentative_index), phis(:, tentative_index),nr,dr,Re,Delta,DTN,Fr,...
                        We,Ma,zs,IntMat(numl(tentative_index)-2,:),angleDropMP,Cang,Dr,RvTent);
                    if abs(errortan(2,tentative_index+1)) < abs(errortan(1,tentative_index+1))
                        numlTent = numl(tentative_index)-1;
                        etaTent = etaprob(:,2);
                        phiTent = phiprob(:,2);
                        psNew = psprob(:,2);
                        zTent = zprob(2);
                        vzTent = vzprob(2);
                    else
                        tvec = [tvec(1:tentative_index),tvec(tentative_index)/2+tvec(tentative_index+1)/2,tvec(tentative_index+1:end)];
                        tentative_index = tentative_index-1; 
                        reduc = 1;
                    end
                else
                    numlTent = numl(tentative_index);
                    etaTent = etaprob(:,3);
                    phiTent = phiprob(:,3);
                    psNew = psprob(:,3);
                    zTent = zprob(3);
                    vzTent = vzprob(3);
                end
            else
                co = find(numl(tentative_index:-1:1)~=numl(tentative_index)-1,1);
                [etaprob(:,2),phiprob(:,2),zprob(2),vzprob(2),psprob(1:numl(tentative_index)-1,2),errortan(2,tentative_index+1)] = ...    
                    solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index)-1,dt,z(tentative_index),vz(tentative_index),etas(:, tentative_index), phis(:, tentative_index),nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(numl(tentative_index)-1,:),angleDropMP,Cang,Dr,RvTent);
                if abs(errortan(2,tentative_index+1)) < 4
                    co = find(numl(tentative_index:-1:1)~=numl(tentative_index)-2,1);
                    [~,~,~,~,~,errortan(1,tentative_index+1)] = ...
                        solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index)-2,dt,z(tentative_index),vz(tentative_index),etas(:, tentative_index), phis(:, tentative_index),nr,dr,Re,Delta,DTN,Fr,...
                        We,Ma,zs,IntMat(numl(tentative_index)-2,:),angleDropMP,Cang,Dr,RvTent);
                    if abs(errortan(2,tentative_index+1)) < abs(errortan(1,tentative_index+1))
                        numlTent = numl(tentative_index)-1;
                        etaTent = etaprob(:,2);
                        phiTent = phiprob(:,2);
                        psNew = psprob(:,2);
                        zTent = zprob(2);
                        vzTent = vzprob(2);
                    else
                        tvec = [tvec(1:tentative_index),tvec(tentative_index)/2+tvec(tentative_index+1)/2,tvec(tentative_index+1:end)];
                        tentative_index = tentative_index-1; 
                        reduc = 1;
                    end
                else
                    tvec = [tvec(1:tentative_index),tvec(tentative_index)/2+tvec(tentative_index+1)/2,tvec(tentative_index+1:end)];
                    tentative_index = tentative_index-1; 
                    reduc = 1;
                end
            end

            if ll == 100 && reduc == 0
                tvec = [tvec(1:tentative_index),tvec(tentative_index)/2+tvec(tentative_index+1)/2,tvec(tentative_index+1:end)];
                tentative_index = tentative_index-1;
                reduc = 1;
            end

            % check whether warning was raised during loop
            [~, warnId] = lastwarn();

            % If there was a warning of dt is too small, step back one step
            if strcmp(warnId, 'MATLAB:nearlySingularMatrix') == true || dt*T_unit < 1e-10
                lastwarn('', '');

                % bc if reduc == 1, tentative index was already lowered
                if reduc == 1; tentative_index = tentative_index + 1; end

                % eliminate time variables that are too close to that time
                tvec(tvec >= tvec(tentative_index) & tvec < (tvec(tentative_index) + dtb)) = [];
                g = @(a, b) a + 0.6 * (b-a);
                tvec = [tvec(1:(tentative_index-1)),...
                        g(tvec(tentative_index-1),tvec(tentative_index)),...
                        tvec(tentative_index:end)];
                tentative_index = tentative_index - 2;
                reduc = 1;

                % reinstating old variables
                %etao = eta_old;
                %phio = phi_old;
                ps_accepted = ps_old;
                previous_conditions{2} = previous_conditions{1};
                previous_conditions{1} = old_conditions;

                going_back = going_back + 1;
                disp("Warning detected. Will proceed with going back"); 
                if going_back > 50
                    %exit = true;
                    PROBLEM_CONSTANTS.problem_flag = true;
                    error('Went back too many times. Stopping execution');
                    
                end
            % verifying convergence of the pressure field
            elseif reduc == 0 %if there was no reduction of time step
                if norm(psNew,1) == 0
                    B_l_ps_new = zeros(1, N);
                else
                    nb_contact_points = numlTent; %nlmaxTent-find(flipud(psNew),1)+1;%Number of nodes in which the pressure needs to be integrated
                    %against each harmonic 
                    if PROBLEM_CONSTANTS.linear_on_theta == true
                        if nb_contact_points > 1
                            contactAngle = (1.5 * thetaVec(nb_contact_points) - 0.5*thetaVec(nb_contact_points-1));
                        else
                            contactAngle = (thetaVec(2) + thetaVec(1))/2;
                        end
                        angles = linspace(contactAngle, ...% + (thetaVec(nb_contact_points) - thetaVec(nb_contact_points-1))/2
                                    pi, (nb_contact_points +1) * PROBLEM_CONSTANTS.interpolation_number);
                        values = r_from_spherical(angles, oscillation_amplitudes(:, tentative_index));
                        f = @(r) interp1(dr*(0:(nb_contact_points)), ...
                            [psNew(1:nb_contact_points)', 0], r, 'linear',  0);
                        values = f(values);
                        B_l_ps_new = custom_project_amplitudes(angles, values, N, NaN, NaN);

                    else
                        % Defining where the pressure distribution will be
                        % projected
                        if nb_contact_points == length(thetaVec)
                            angles = [thetaVec(1:(nb_contact_points)), (2 * thetaVec(nb_contact_points) - thetaVec(nb_contact_points-1))];
                        else
                            angles = thetaVec(1:(nb_contact_points+1));
                        end
                        f = @(thetas) interp1(angles, [psNew(1:nb_contact_points)', 0], thetas, 'linear',  0); 
                        endpoints = [angles(end), angles(1)];
                        B_l_ps_new = project_amplitudes(f, N, endpoints, PROBLEM_CONSTANTS, true); 
                    end
                end

                [amplitudes_new, velocities_new] = solve_ODE_unkown(nan, B_l_ps_new, dt, ...
                    previous_conditions, PROBLEM_CONSTANTS);

                errorP = norm(amplitudes_tent-amplitudes_new)/norm(amplitudes_tent);

                if PROBLEM_CONSTANTS.DEBUG_FLAG == true
                    fprintf("Inside ll: %0.2g, errP: %0.5g \n", ll, errorP);
                    disp(errortan(:, tentative_index+1)');
                end

                if errorP < tolP % Finally accept solution

                    % Saving old variables in case we need to go back
                    %numlOld = numl(tentative_index);
                    %eta_old = eta_accepted;
                    %phi_old = phi_accepted;
                    ps_old = ps_accepted;
                    %z_old = z(tentative_index);
                    %vz_old = vz(tentative_index);
                    %amplitudes_old = oscillation_amplitudes(:, tentative_index);
                    %B_l_ps_old = pressure_amplitudes(:, tentative_index);
                    %Rv_old = Rv(tentative_index);
                    old_conditions = previous_conditions{1};

                    % Storing accepted variables
                    numl(tentative_index+1) = numlTent;
                    eta_accepted = etaTent;
                    phi_accepted = phiTent;
                    ps_accepted = psNew;
                    z(tentative_index+1) = zTent;
                    vz(tentative_index+1) = vzTent;
                    %#---
                    oscillation_amplitudes(:, tentative_index + 1) = amplitudes_new;
                    pressure_amplitudes(:, tentative_index + 1)    = B_l_ps_new;
                    Rv(tentative_index+1) = zs_from_spherical(pi, amplitudes_new);
                    %amplitudes_velocities_old = velocities_new;


                    previous_conditions{1} = previous_conditions{2};
                    previous_conditions{2} = struct("deformation_amplitudes", amplitudes_new, ...
                        "deformation_velocities", velocities_new, ...
                        "dt", dt, "nb_harmonics", N, "pressure_amplitudes", B_l_ps_new, ...
                        "current_time", previous_conditions{1}.current_time + dt, ...
                        "center_of_mass", z(tentative_index+1), "center_of_mass_velocity", vz(tentative_index+1), ...
                        "nb_contact_points", numlTent);

                    nlmax(tentative_index+1) = nlmaxTent;
                    etaOri(tentative_index+1) = eta_accepted(1);
                    etas(:, tentative_index + 1) = eta_accepted;
                    phis(:, tentative_index + 1) = phi_accepted;

                    if t >= tvecOri(current_to_save) || t < 10*dtb
                        %etaMatPer(:,current_to_save) = eta_accepted;
                        %phiMatPer(:,current_to_save) = phi_accepted;
                        %psMatPer{current_to_save} = ps_accepted;
                        indexes_to_save(current_to_save) = tentative_index;
                        current_to_save = current_to_save + 1;
                    end

                    if (zTent > 1.5 && numlTent == 0) || (vz(tentative_index+1) < 0 && vz(tentative_index) >=0  && t > 50 * dtb)
                        tend = t;
                    end

                    %etao = eta_accepted;
                    %phio = phi_accepted;
                    %pso = ps_accepted;

                    if PROBLEM_CONSTANTS.DEBUG_FLAG == true
                        xs = dr*(0:nlmax(tentative_index+1)-1);
                        zsplot = zs(1:nlmax(tentative_index+1))+RvTent+z(tentative_index+1);
                        plot([-fliplr(xs(2:end)),xs],[flipud(zsplot(2:end));zsplot],'k','Linewidth',2);
                        hold on
                        thetaplot = linspace(0, thetaVec(end), 200);%-%-0:thetaVec(end)/400:thetaVec(end);
                        %-%-xsTop = xsoftheta(thetaplot,A2New,A3New);
                        %-%-zsTop = zsoftheta(thetaplot,A2New,A3New);
                        zsTop = zs_from_spherical(thetaplot, amplitudes_new);
                        xsTop = r_from_spherical(thetaplot, amplitudes_new); 
                        plot([-xsTop(end:-1:2), xsTop],[zsTop(end:-1:2), zsTop]+zTent,'k','Linewidth',2);
                        width = min(nr, 200);
                        plot([-fliplr(xplot(2:width)),xplot(1:width)],[flipud(eta_accepted(2:width));eta_accepted(1:width)],'LineWidth',2);
                        hold off
                        axis equal
                        title(sprintf('   t = %0.3f, nl = %d', t, numl(tentative_index+1)),'FontSize',16);
                        grid on
                        %set(gca,'xlim',[-6 6])
                        drawnow;
                    end
                else %i.e. didn't attain convergence yet

                    
                    %psTent = psNew;

                    %#---
                    amplitudes_tent = amplitudes_new;
                    %B_l_ps_tent = B_l_ps_new;
                    %#---

                    RmaxTent = r_from_spherical(maximum_contact_radius(amplitudes_tent), amplitudes_tent);

                    %RmaxTent = r_from_spherical(theta_max_radius(amplitudes_tent), amplitudes_tent);
                    %-%-RmaxTent = xsoftheta(mod(abs(thetaMax(A2Tent,A3Tent)),pi),A2Tent,A3Tent);
                    nlmaxTent = floor(RmaxTent/dr)+1;

                    thetaVec = theta_from_cylindrical(dr*(0:(nlmaxTent-1)), oscillation_amplitudes(:, tentative_index));

                    %-%-RvTent = zsoftheta(pi,A2Tent,A3Tent);

                    %#---
                    RvTent = zs_from_spherical(pi, amplitudes_tent);
                    zs(1:nlmaxTent) = zs_from_spherical(thetaVec, amplitudes_tent)' - RvTent;
                    zs((nlmaxTent+1):nr) = Inf;
                    %#---

                    %-%-zs(1:nlmaxTent) = zsoftheta(thetaVec,A2Tent,A3Tent)-RvTent;
                    %-%-zs(nlmaxTent+1:end) = Inf;

                    %#---
                    tanDrop = calculate_tan( dr * (1:nlmaxTent) - dr/2, amplitudes_tent)';
                    angleDropMP(1:(nlmaxTent)) = atan(tanDrop(1:(nlmaxTent))); %% TODO: Check these for boundary points
                    %finding angle at tangent direction

                    psprob = zeros(nlmaxTent,5);%zeroing the vector of potential pressures
                end

            end %it time step was not reduced
        end % inner while   
    end % Outer while




    
    for ii = 1:length(savingvarNames)
       variableValues{ii} = eval(savingvarNames{ii}); 
    end

    indexes_to_save = indexes_to_save(1:(current_to_save-1));
    results_saver("", indexes_to_save, variableValues, savingvarNames);

catch ME
    for ii = 1:length(savingvarNames)
       variableValues{ii} = eval(savingvarNames{ii}); 
    end

    indexes_to_save = indexes_to_save(1:(current_to_save-1));
    results_saver("errored_", indexes_to_save, variableValues, savingvarNames);
    
    
        
    fprintf("Couldn't run simulation with the following parameters: \n Velocity: %g \n Modes: %g \n", ...
        U0, N); 
    a = datetime('now'); a.Format = 'yyyyMMddmmss';
    save(sprintf("error_logU0=%g-%s.mat", U0, a),'ME');
end % end while catch

simul_time = toc(tstart);
%simul_time = simul_time - tstart;

save('ProblemConditions.mat', "T", "N", "U0", "Ang", "Re", "Fr", "We", ...
"WeS", "Cang", "tend", "nsteps", "dtb", "L_unit", "T_unit", "M_unit", ...
"PROBLEM_CONSTANTS", "simul_time");
mypwd = split(pwd, "1_code"); mypwd = mypwd{2};
fprintf("Finished simulation on %s. Time elapsed: %0.2f minutes\n", mypwd, simul_time/60);
cd(currfold)


end

function results_saver(prefix, indexes, variables, variableNames)
    if indexes(2) == 1
        indexes = indexes(2:end); 
    end
    for ii = 1:length(variables)
       var = variables{ii};
       switch length(size(var))
           case 1
               var = var(indexes);
           case 2
               if iscell(var)
                   var = var{:, indexes};
               else
                   var = var(:, indexes);
               end
       end
       stru = struct(variableNames{ii}, var);
       save(sprintf('%s%s.mat', prefix, variableNames{ii}), '-struct', 'stru');
    end
    
end

function out = getVarName(var)
    out = inputname(1);
end
