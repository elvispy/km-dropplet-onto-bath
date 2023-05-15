%clear
close all
%clc

tstart = tic;
tic
%data in cgs
tmax = 100;

%load('runNumber.mat','runNumber'); 
runNumber = 0; %-%-

if exist('z.mat', 'file') == 2
   % error("Exporting data is going to be overwritten. Please re-allocate files to avoid loss of data");
end

if runNumber == 0
    U0 = 38; save('U0.mat','U0')%impact velocity in cm/s (unit of velocity for the problem)
    Ang = 180; save('Ang.mat','Ang') %contact angle to be imposed
   
    cd ..
    load('Ro.mat','Ro')%Sphere's radius in CGS
    
    cd ..
    load('rhoS.mat','rhoS')%Sphere density
    load('sigmaS.mat')%Sphere's surface tension
    
    cd ..
    load('rho.mat','rho')
    load('sigma.mat','sigma')
    load('nu.mat','nu')
    load('muair.mat')
    load('g.mat','g') %gravitational constant
    
    cd ..
    load('D.mat')%Domain diameter in units of droplet radii
    load('quant.mat')%number of dr's contained in an undeformed dropelt radius
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
    
    cd(['ImpDefCornerAng',num2str(Ang),'U',num2str(U0)])

    tiempoComp = zeros(1,10); %just to check how long it takes to solve the first ten saving intervals
    
    % #--- 
    N = 40; % Number of harmonics contributing to the oscillation
    % #---0
    
    %Unit of time
    T = Ro/U0; save('T.mat','T')%base time is seconds

    %Dimensionless numbers that depend on U0
    Re = Ro*U0/nu; save('Re.mat','Re')
    Fr = U0^2/(g*Ro); save('Fr.mat','Fr')
    We = rho*Ro*U0^2/sigma; save('We.mat','We')
    WeSB = rhoS*Ro*U0^2/sigma;save('WeSB.mat','WeSB')
    WeS  = rhoS*Ro*U0^2/sigmaS;save('WeS.mat','WeS') %This name may not be the best, the surface tension is that of the 
    %bath at least in one place
    Cang = (Ang/180)*pi; save('Cang.mat','Cang')%contact angle to be imposed
    
    %Physical parameters
    tend = 9; save('tend.mat','tend')%Earliest possible end of simulation in characteristic units
    
    %Inintial conditions for the fluid
    t = 0;
    etao = zeros(nr,1); %initial surface elevation
    phio = zeros(nr,1); %initial surface potential

    %Numerical Simulation parameters
    nsteps = 1000; save('nsteps.mat','nsteps')%minimum number of timesteps in one unit of time
    dtb = 1/nsteps; save('dtb.mat','dtb')%basic timestep (gets halved as needed over impacts)
    steps = ceil((tend-t)/dtb); %estimated minimum number of timesteps
    
    %Zeroing result storing variables
    etaOri = zeros(1,steps+1);%height of the surface below the south pole
    z = zeros(1,steps+1);%height of the centre of mass
    vz = zeros(1,steps+1);%speed of the centre of mass
    numl = zeros(1,steps+1);%number of pressed mesh points at each time step
    tvec = t:(dtb):tend+1; %vector of times assuming no refinement has happened
    %plus some extra time just in case the simulation needs to run longer
    % #--- 
    oscillation_amplitudes = zeros(N, steps + 1); % Variable to store
    Rv = -ones(1, steps+1);
    % the time dependent amplitude of all the SH
    oscillation_velocities = zeros(N, steps+1);
    % #--- 
    %-%-A2 = zeros(1,steps+1);%Variable to store the time dependent amplitude of the 2nd soherical harmonic (SH) mode
    %-%-V2 = zeros(1,steps+1);%Variable to store the time dependent velocity  of the 2nd SH mode
    %-%-A3 = zeros(1,steps+1);%Variable to store the time dependent amplitude of the 3nd SH mode
    %-%-V3 = zeros(1,steps+1);%Variable to store the time dependent velocity  of the 3nd SH mode
    nlmax = zeros(1,steps+1);%Variable to store the number of nodes spanned by the deformed droplet
    
    tolP = 1E-6; save('tolP.mat','tolP')%error tolerance for the pressure field and deformation 
    
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
    vz(1) = -1; %Initial velocity of the CoM in dimesionless units
    dt = tvec(2) - tvec(1);
    
    current_conditions = struct("deformation_amplitudes", amplitudes_old, ...
        "deformation_velocities", amplitudes_velocities_old, ...
        "pressure_amplitudes", B_l_ps_old, "dt", dt, "nb_harmonics", N,  ...
        "current_time", 0, ...
        "center_of_mass", z(1), "center_of_mass_velocity", vz(1), ...
        "nb_contact_points", 0);
    
    previous_conditions = {current_conditions, current_conditions}; 
    
    f = @(n)  sqrt(n .* (n+2) .* (n-1) / WeS);
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
    
    jj = 0;%iteration counter

    errortan = zeros(5,steps+1);%tangency error recorder

    %intial guess for the surface profile at the next time step
%     eta1 = etao;%these are here just to satisfy a silly thing in the first loop
%     phi1 = phio;
    ps1 = [];
    
    % #---
    
    %If there were some initial pressure acting on the surface and sphere I
    %would have to change this bit here to reflect the presure distribution
    
    %zeroing variable that records each part of the sequence of surface states
    etaMatPer = zeros(length(etao),nsteps); % ?? Why nsteps
    phiMatPer = zeros(length(phio),nsteps);
    psMatPer = cell(1,nsteps); % ??? why cell
    %Storing initial surface state
    etaMatPer(:,1) = etao;
    phiMatPer(:,1) = phio;
    psMatPer{1} = zeros(quant+1,1);
    
    %zeroing the ceiling functions
    zs = zeros(nr,1);
    
    jj1 = 1; %partial results savings  counter
    
    PROBLEM_CONSTANTS = struct("froude_nb", Fr, "weber_nb", WeS, ...
        "nb_harmonics", N, ...
        "omegas_frequencies", omegas_frequencies, ...
        "spatial_tol", dr, ...
        "DEBUG_FLAG", true, ...
        "Ra", Ra);
                                %"pressure_unit", pressure_unit, ...
                                %"CM", 9, ...
                                %"PG", 2, ...
                                %"KILL_OUTSIDE", true, ...
                                %"wigner3j", {precomputed_wigner(harmonics_qtt)}, ...
        
end
%% hey
% 
% elseif runNumber>0
%     load('U0.mat','U0')%impact velocity in cm/s (unit of velocity for the problem)
%     load('Ang.mat','Ang') %contact angle to be imposed
%     
%     cd ..
%     load('Ro.mat','Ro')%Sphere's radius in CGS
%     load('D.mat')%Domain diameter in units of droplet radii
%     load('quant.mat')%number of dr's contained in an undeformed dropelt radius
%     
%     cd ..
%     load('rhoS.mat','rhoS')%Sphere density
%     load('sigmaS.mat')%Sphere's surface tension
%     
%     cd ..
%     load('sigma.mat','sigma')
%     
%     %loading fluid properties
%     load('rho.mat','rho')
%     load('nu.mat','nu')
%     load('g.mat','g') %gravitational constant
%     
%     cd(['RhoS',num2str(rhoS*1000),'SigmaSoverSigma',num2str(sigmaS/sigma)])
%     load('Ma.mat','Ma')%Dimensionless mass of sphere
%     load('Ra.mat','Ra')%Density ratio
%     
%     cd(['R0',num2str(Ro*10000),'mmD',num2str(D),'xRoQuant',num2str(quant)])
%     load('dropmass.mat','dropmass')%Mass of droplet in cgs
%     load('nr.mat','nr')
%     load('dr.mat','dr')
%     load('Delta.mat','Delta')
%     load('IntMat.mat','IntMat')
%     load('DTNnew345nr2500D100refp10.mat','DTNnew345')
%     DTN = DTNnew345;
%     clear DTNnew345
%     load('xplot.mat','xplot');%I might remove or relocate this
%     
%     cd(['ImpDefCornerAng',num2str(Ang),'U',num2str(U0)])
%     
%     %Unit of time
%     load('T.mat','T')%base time is seconds
% 
%     %Dimensionless numbers that depend on U0
%     load('Re.mat','Re')
%     load('Fr.mat','Fr')
%     load('We.mat','We')
%     load('WeSB.mat','WeSB')
%     load('WeS.mat','WeS') %This name may not be the best, the surface tension is that of the 
%     %bath at least in one place
%     load('Cang.mat','Cang')%contact angle to be imposed
%     
%     %Physical parameters
%     load('tend.mat','tend')%Earliest possible end of simulation in characteristic units
%     
%     %Inintial conditions for the fluid
%     t = 0;
%     etao = zeros(nr,1); %initial surface elevation
%     phio = zeros(nr,1); %initial surface potential
% 
%     %Numerical Simulation parameters
%     load('nsteps.mat','nsteps')%minimum number of timesteps in one unit of time
%     load('dtb.mat','dtb')%basic timestep (gets halved as needed over impacts)
%     
%     %Storing variables
%     load('etaOri.mat','etaOri')%height of the surface below the south pole
%     load('z.mat','z') %height of the centre of mass
%     load('vz.mat','vz')%speed of the centre of mass
%     load('numl.mat','numl')%number of pressed mesh points at each time step
%     load('tvec.mat','tvec')
%     load('A2.mat','A2')%Variable to store the time dependent amplitude of the 2nd soherical harmonic (SH) mode
%     load('V2.mat','V2')%Variable to store the time dependent velocity  of the 2nd SH mode
%     load('A3.mat','A3')%Variable to store the time dependent amplitude of the 3nd SH mode
%     load('V3.mat','V3')%Variable to store the time dependent velocity  of the 3nd SH mode
%     load('nlmax.mat','nlmax');%Variable to store the number of nodes spanned by the deformed droplet
%     
%     load('tolP.mat','tolP')%error tolerance for the pressure field and deformation
%     
%     %Drop oscillation frequencies
%     load('omega2.mat','omega2')%Angular frequency of 2nd SH mode
%     load('omega3.mat','omega3')%Angular frequency of 3rd SH mode
%     load('P2.mat','P2')%Matrix whose columns are the eigenvectors of the matrix of 
%     %coefficients of the 2x2 ODE system for the 2nd SH mode
%     load('P2inv.mat','P2inv')%its inverse
%     load('P3.mat','P3')%Matrix whose columns are the eigenvectors of the matrix of 
%     %coefficients of the 2x2 ODE system for the 3rd SH mode
%     load('P3inv.mat','P3inv')%its inverse
%     
%     %Inintial conditions
%     load(['tstop',num2str(runNumber),'.mat'],'tstop') 
%     t = tstop;
%     load(['jjstop',num2str(runNumber),'.mat'],'jjstop')
%     jj = jjstop;
%     
%     %Initial conditions for the drop (shape and deformation speed for each SH mode)
%     A2old = A2(jj);
%     V2old = V2(jj);
%     A3old = A3(jj);
%     V3old = V3(jj);
%     Y2Tent = zeros(2,1);%Intialising canonical variables for drop oscillation integration
%     Y3Tent = zeros(2,1);
%     Y2New = zeros(2,1);%Intialising canonical variables for drop oscillation integration
%     Y3New = zeros(2,1);
% 
%     load('errortan.mat','errortan')%tangency error recorder
% 
%     %Initial surface state
%     load(['etao',num2str(runNumber),'.mat'],'etao')
%     load(['phio',num2str(runNumber),'.mat'],'phio')
%     load(['pso',num2str(runNumber),'.mat'],'pso')
%     ps1 = pso;
%     
%     %I need to modify this part to account for a non-zero
%     %pso%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%
%     B2old = 0;%2nd SH component of initial pressure
%     B3old = 0;%3rd SH component of initial pressure
%     
%     %zeroing variable that records each part of the sequence of surface states
%     etaMatPer = zeros(length(etao),nsteps);
%     phiMatPer = zeros(length(phio),nsteps);
%     psMatPer = cell(1,nsteps);
%     %Storing initial surface state
%     etaMatPer(:,1) = etao;
%     phiMatPer(:,1) = phio;
%     psMatPer{1} = pso;
%     
%     %zeroing the ceiling functions
%     zs = zeros(nr,1);
%     % I need to review this zs part
%     
%     jj1 = 1; %partial results savings  counter
% 
% elseif runNumber < 0
%     exit
% end

%% Main Loop
while (t<tend) %#-- || jj1>.5) 
%     if toc > 120
%         break;
%     end
    jj = jj+1;
    t = tvec(jj+1);
    dt = t - tvec(jj);
    %zeroing the tentative solution variables
    etaprob = zeros(nr,5);
    phiprob = zeros(nr,5);
    vzprob  = zeros(1,5);
    zprob   = zeros(1,5);
    errortan(:,jj+1) = 4*ones(5,1);
    
    psTent = ps1; %Tentative pressure distribution (we start with the previous pressure)
    
    %Shape at the start of the time step
    %x that corresponds to the max pressed radius
    %-%-RmaxOld = xsoftheta(mod(abs(thetaMax(A2(jj),A3(jj))),pi),A2(jj),A3(jj)); 
    %#---
    %function_amplitudes = oscillation_handle(oscillation_amplitudes(:, jj));
    %function_amplitudes_prime = oscillation_handle_prime(oscillation_amplitudes(:, jj));
    %function_amplitudes_prime_prime = oscillation_handle_prime_prime(oscillation_amplitudes(:, jj));
    
    RmaxOld = r_from_spherical(maximum_contact_radius(oscillation_amplitudes(:, jj)), oscillation_amplitudes(:, jj));
  
    %(i.e. where the tangent plane to the droplet is vertical)
    nlmax(jj) = floor(RmaxOld/dr)+1;%max number of contact points
    % #---
    %nlmax_dummy = floor(max_radius/dr)+1;
    
    thetaVec = theta_from_cylindrical(dr*(0:(nlmax(jj)-1)), oscillation_amplitudes(:, jj)); % zeros(1,nlmax(jj));%initialising vector of angles of pressed positions
%     thetaVec(1) = pi;%bottom of droplet
%     for ii = 2:nlmax(jj)
%         %-%-thetaVec(ii) = mod(abs(thetaofxs(dr*(ii-1),A2(jj),A3(jj),thetaVec(ii-1))),pi); % #--- I think we could make this analytic.
%         thetaVec(ii) = theta_from_cylindrical(dr*(ii-1), oscillation_handle(oscillation_amplitudes(:, jj)), ...
%                         oscillation_handle_prime(oscillation_amplitudes(:, jj)), thetaVec(ii-1));
%     end
    % xi = cos(thetaVec); %Variable change to use the Legendre polynomials
    %finding tentative B2 and B3, i.e. projection of the pressure onto the
    %Spherical Harmonic modes
    %downward incline at origin
    if norm(psTent,1) == 0
        %B2Tent = 0;
        %B3Tent = 0;
        B_l_ps_tent = zeros(1, N);
    else
        nb_contact_points = nlmax(jj)-find(flipud(psTent),1)+1; %Number of nodes contact points%
        %needs to be integrated against SH modes
        
        %%#---
        % We have that B_l = (2l+1)/2 * int(Ps * P_m * sin(theta))
%         B_l_ps_tent = arrayfun(@(m) find_harmonic_coefficient(thetaVec(1:(nb_contact_points+1)), ...
%                     psTent(1:nb_contact_points), m, nb_contact_points, LEGENDRE_POLYNOMIALS{m}), 1:N);
        
        f = @(thetas) interp1(thetaVec(1:(nb_contact_points+1)), [psTent(1:nb_contact_points)', 0], thetas, 'linear',  0); 
        endpoints = [thetaVec(nb_contact_points+1), thetaVec(1)]; %TODO: Check integration ends
        B_l_ps_tent = project_amplitudes(f, N, endpoints, PROBLEM_CONSTANTS, true);

        
    end
    
    
    %#--- Solving the ODE
    %C_tent = zeros(2, N);
%     [amplitudes_tent, ~] = solve_EDO(N, dt, Ra, ...
%         omegas_frequencies, B_l_ps_old, B_l_ps_tent, ODE_inverse_matrices, ...
%         ODE_matrices, amplitudes_old, amplitudes_velocities_old);%zeros(1, N);
    [amplitudes_tent, velocities_tent] = solve_ODE_unkown(nan, B_l_ps_tent, dt, ...
        previous_conditions, PROBLEM_CONSTANTS);
    %#---
    
    
    %#---
%     function_amplitudes = oscillation_handle(amplitudes_tent);
%     function_amplitudes_prime = oscillation_handle_prime(amplitudes_tent);
%     function_amplitudes_prime_prime = oscillation_handle_prime_prime(amplitudes_tent);
    
%     RmaxTent = r_from_spherical(theta_max_radius(function_amplitudes, function_amplitudes_prime, ...
%         function_amplitudes_prime_prime), function_amplitudes);
    RmaxTent = r_from_spherical(maximum_contact_radius(oscillation_amplitudes(:, jj)), oscillation_amplitudes(:, jj));
  
    %RmaxTent = r_from_spherical(theta_max_radius(amplitudes_tent), amplitudes_tent);
    %Tentative shape
    %-%-RmaxTent = xsoftheta(mod(abs(thetaMax(A2Tent,A3Tent)),pi),A2Tent,A3Tent);
    nlmaxTent = floor(RmaxTent/dr)+1;
    thetaVec  = theta_from_cylindrical(dr*(0:(nlmaxTent-1)), oscillation_amplitudes(:, jj)); % zeros(1,nlmaxTent);
%     thetaVec(1) = pi;
%     for ii = 2:nlmaxTent
%         %-%-thetaVec(ii) = mod(abs(thetaofxs(dr*(ii-1),A2Tent,A3Tent,thetaVec(ii-1))),pi);
%         thetaVec(ii) = theta_from_cylindrical(dr*(ii-1), oscillation_handle(amplitudes_tent), ...
%                         oscillation_handle_prime(amplitudes_tent), thetaVec(ii-1));
%     end
%     xi = cos(thetaVec);
%     
    %-%-RvTent = zsoftheta(pi,A2Tent,A3Tent);%height of the south pole with respect to the centre of mass
    %#---
    RvTent = zs_from_spherical(pi, amplitudes_tent);
    zs(1:nlmaxTent) = zs_from_spherical(thetaVec, amplitudes_tent)' - RvTent; %TODO: Check that matrix dimensions agree.
    zs((nlmaxTent+1):nr) = Inf;
    %#---
    %-%-zs(1:nlmaxTent) = zsoftheta(thetaVec,A2Tent,A3Tent)-RvTent;%Height of the bottom boundary of the
    %droplet with respect to the south pole
    %-%-zs(nlmaxTent+1:nr) = Inf;%height outside of the droplet shadow
    
    %finding angle at tangent direction %may be improved by taking tangent
    %at midle point
    
    %#---
    tanDrop = calculate_tan( dr * (1:nlmaxTent) - dr/2, amplitudes_tent)';
    angleDropMP(1:(nlmaxTent)) = atan(tanDrop(1:(nlmaxTent)));
    
    %for ii  = 2:nlmaxTent
    %    tanDrop_dummy(ii) = calculate_tan(dr* (ii-1/2), amplitudes_tent, legendrepDerivatives);
    %end
    
%     thetaVecMP = zeros(nlmaxTent,1);
%     thetaVecMP(1) = mod(abs(thetaofxs(.5*dr,A2Tent,A3Tent,pi)),pi);
%     for ii = 2:nlmaxTent-1
%         thetaVecMP(ii) = mod(abs(thetaofxs(dr*(ii-.5),A2Tent,A3Tent,thetaVecMP(ii-1))),pi);
%     end
%     thetaVecMP(nlmaxTent) = mod(abs(thetaofxs(dr*(nlmaxTent-1),A2Tent,A3Tent,thetaVecMP(nlmaxTent-1))),pi);
%     ysp = ysprimeoftheta(thetaVecMP,A2Tent,A3Tent);%derivative of y with respect to theta at nodes
%     xsp = xsprimeoftheta(thetaVecMP,A2Tent,A3Tent);%derivative of x with respect to theta at the nodes
%     angleDropMP = zeros(nlmaxTent,1);
%     if xsp(end) == 0 %if last point's tangent is infinity
%         tanDrop(1:nlmaxTent-1) = ysp(1:end-1)./xsp(1:end-1);
%         angleDropMP(1:nlmaxTent-1) = atan(tanDrop(1:nlmaxTent-1));
%         angleDropMP(nlmaxTent) = pi/2;
%     else
%         tanDrop(1:nlmaxTent) = ysp./xsp;
%         angleDropMP(1:nlmaxTent) = atan(tanDrop(1:nlmaxTent));
%     end
    
    psprob = zeros(nlmaxTent,5);%zeroing the vector of potential pressures
    errorP = 1; %error in the pressure field and amplitude of modes
    reduc = 0; %indicator of whether there was a reduction in the time-step size or not
    ll = 0; % Limiting while loop to 100 iterations
    while abs(errorP)>=tolP && reduc == 0 
        ll = ll + 1;
        
        if numl(jj) < .5 %i.e. if previously in flight (I need to define this as integer)
            [etaprob(:,3),phiprob(:,3),zprob(3),vzprob(3),errortan(3,jj+1)] = ...
                solveDD0(dt,z(jj),vz(jj),etao,phio,nr,Re,Delta,DTN,Fr,We,zs,RvTent);
            if abs(errortan(3,jj+1))<.5
                numlTent = 0;
                etaTent = etaprob(:,3);
                phiTent = phiprob(:,3);
                psNew = zeros(nlmaxTent,1);
                zTent = zprob(3);
                vzTent = vzprob(3);
            else
                co = find(numl(jj:-1:1)~=1,1);
                [etaprob(:,4),phiprob(:,4),zprob(4),vzprob(4),psprob(1,4),errortan(4,jj+1)] = ...
                    solvenDDCusp(numl(jj-co+1),1,dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(1,:),angleDropMP,Cang,WeSB,RvTent);
                co = find(numl(jj:-1:1)~=2,1);
                [~,~,~,~,~,errortan(5,jj+1)] = ...    
                    solvenDDCusp(numl(jj-co+1),2,dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(2,:),angleDropMP,Cang,WeSB,RvTent);
                if abs(errortan(4,jj+1)) < abs(errortan(5,jj+1))
                    numlTent = 1;
                    etaTent = etaprob(:,4);
                    phiTent = phiprob(:,4);
                    psNew = psprob(:,4);
                    zTent = zprob(4);
                    vzTent = vzprob(4);
                else
                    tvec = [tvec(1:jj),tvec(jj)/2+tvec(jj+1)/2,tvec(jj+1:end)];
                    jj = jj-1;
                    reduc = 1;
                end
            end
        elseif numl(jj)>.5 && numl(jj)<1.5 % i.e. the last number of contact points was 1
            [etaprob(:,2),phiprob(:,2),zprob(2),vzprob(2),errortan(2,jj+1)] = ...
                solveDD0(dt,z(jj),vz(jj),etao,phio,nr,Re,Delta,DTN,Fr,We,zs,RvTent);
            if abs(errortan(2,jj+1))<.5
                numlTent = 0;
                etaTent = etaprob(:,2);
                phiTent = phiprob(:,2);
                psNew = zeros(nlmaxTent,1);
                zTent = zprob(2);
                vzTent = vzprob(2);
            else
                co = find(numl(jj:-1:1)~=1,1);
                [etaprob(:,3),phiprob(:,3),zprob(3),vzprob(3),psprob(1,3),errortan(3,jj+1)] = ...                     
                    solvenDDCusp(numl(jj-co+1),1,dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(1,:),angleDropMP,Cang,WeSB,RvTent);
                co = find(numl(jj:-1:1)~=2,1);
                [etaprob(:,4),phiprob(:,4),zprob(4),vzprob(4),psprob(1:2,4),errortan(4,jj+1)] = ...    
                    solvenDDCusp(numl(jj-co+1),2,dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(2,:),angleDropMP,Cang,WeSB,RvTent);
                if abs(errortan(3,jj+1)) < abs(errortan(4,jj+1))
                    numlTent = 1;
                    etaTent = etaprob(:,3);
                    phiTent = phiprob(:,3);
                    psNew = psprob(:,3);
                    zTent = zprob(3);
                    vzTent = vzprob(3);
                else
                    co = find(numl(jj:-1:1)~=3,1);
                    [~,~,~,~,~,errortan(5,jj+1)] = ...
                        solvenDDCusp(numl(jj-co+1),3,dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                        We,Ma,zs,IntMat(3,:),angleDropMP,Cang,WeSB,RvTent);
                    if abs(errortan(4,jj+1)) < abs(errortan(5,jj+1))
                        numlTent = 2;
                        etaTent = etaprob(:,4);
                        phiTent = phiprob(:,4);
                        psNew = psprob(:,4);
                        zTent = zprob(4);
                        vzTent = vzprob(4);
                    else
                        tvec = [tvec(1:jj),tvec(jj)/2+tvec(jj+1)/2,tvec(jj+1:end)];
                        jj = jj-1; 
                        reduc = 1;
                    end
                end
            end
        elseif numl(jj) > 1.5 && numl(jj) < 2.5 %i.e. the last contact had two points
            [~,~,~,~,errortan(1,jj+1)] = ...
                solveDD0(dt,z(jj),vz(jj),etao,phio,nr,Re,Delta,DTN,Fr,We,zs,RvTent);
            if abs(errortan(1,jj+1))<.5
                tvec = [tvec(1:jj),tvec(jj)/2+tvec(jj+1)/2,tvec(jj+1:end)];
                jj = jj-1; 
                reduc = 1;
            else
                co = find(numl(jj:-1:1)~=2,1);
                [etaprob(:,3),phiprob(:,3),zprob(3),vzprob(3),psprob(1:2,3),errortan(3,jj+1)] = ...    
                    solvenDDCusp(numl(jj-co+1),2,dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(2,:),angleDropMP,Cang,WeSB,RvTent);
                co = find(numl(jj:-1:1)~=1,1);
                [etaprob(:,2),phiprob(:,2),zprob(2),vzprob(2),psprob(1,2),errortan(2,jj+1)] = ...    
                    solvenDDCusp(numl(jj-co+1),1,dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(1,:),angleDropMP,Cang,WeSB,RvTent);
                if abs(errortan(2,jj+1)) < abs(errortan(3,jj+1))
                    numlTent = 1;
                    etaTent = etaprob(:,2);
                    phiTent = phiprob(:,2);
                    psNew = psprob(:,2);
                    zTent = zprob(2);
                    vzTent = vzprob(2);
                else
                    co = find(numl(jj:-1:1)~=3,1);
                    [etaprob(:,4),phiprob(:,4),zprob(4),vzprob(4),psprob(1:3,4),errortan(4,jj+1)] = ...    
                        solvenDDCusp(numl(jj-co+1),3,dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                        We,Ma,zs,IntMat(3,:),angleDropMP,Cang,WeSB,RvTent);
                    if abs(errortan(3,jj+1)) < abs(errortan(4,jj+1))
                        numlTent = 2;
                        etaTent = etaprob(:,3);
                        phiTent = phiprob(:,3);
                        psNew = psprob(:,3);
                        zTent = zprob(3);
                        vzTent = vzprob(3);
                    else
                        co = find(numl(jj:-1:1)~=4,1);
                        [~,~,~,~,~,errortan(5,jj+1)] = ...    
                            solvenDDCusp(numl(jj-co+1),4,dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                            We,Ma,zs,IntMat(4,:),angleDropMP,Cang,WeSB,RvTent);
                        if abs(errortan(4,jj+1)) < abs(errortan(5,jj+1))
                            numlTent = 3;
                            etaTent = etaprob(:,4);
                            phiTent = phiprob(:,4);
                            psNew = psprob(:,4);
                            zTent = zprob(4);
                            vzTent = vzprob(4);
                        else
                            tvec = [tvec(1:jj),tvec(jj)/2+tvec(jj+1)/2,tvec(jj+1:end)];
                            jj = jj-1; 
                            reduc = 1;
                        end
                    end
                end
            end
        elseif numl(jj)>2.5 && numl(jj)<nlmaxTent-1.5 %if the last number of contact points was far from the boundaries
            co = find(numl(jj:-1:1)~=numl(jj),1);
            [etaprob(:,3),phiprob(:,3),zprob(3),vzprob(3),psprob(1:numl(jj),3),errortan(3,jj+1)] = ...    
                solvenDDCusp(numl(jj-co+1),numl(jj),dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                We,Ma,zs,IntMat(numl(jj),:),angleDropMP,Cang,WeSB,RvTent);
            co = find(numl(jj:-1:1)~=numl(jj)-1,1);
            [etaprob(:,2),phiprob(:,2),zprob(2),vzprob(2),psprob(1:numl(jj)-1,2),errortan(2,jj+1)] = ...    
                solvenDDCusp(numl(jj-co+1),numl(jj)-1,dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                We,Ma,zs,IntMat(numl(jj)-1,:),angleDropMP,Cang,WeSB,RvTent);
            if abs(errortan(2,jj+1)) < abs(errortan(3,jj+1))
                co = find(numl(jj:-1:1)~=numl(jj)-2,1);
                [~,~,~,~,~,errortan(1,jj+1)] = ...
                    solvenDDCusp(numl(jj-co+1),numl(jj)-2,dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(numl(jj)-2,:),angleDropMP,Cang,WeSB,RvTent);
                if abs(errortan(2,jj+1)) < abs(errortan(1,jj+1))
                    numlTent = numl(jj)-1;
                    etaTent = etaprob(:,2);
                    phiTent = phiprob(:,2);
                    psNew = psprob(:,2);
                    zTent = zprob(2);
                    vzTent = vzprob(2);
                else
                    tvec = [tvec(1:jj),tvec(jj)/2+tvec(jj+1)/2,tvec(jj+1:end)];
                    jj = jj-1; 
                    reduc = 1;
                end
            else
                co = find(numl(jj:-1:1)~=numl(jj)+1,1);
                [etaprob(:,4),phiprob(:,4),zprob(4),vzprob(4),psprob(1:numl(jj)+1,4),errortan(4,jj+1)] = ...    
                    solvenDDCusp(numl(jj-co+1),numl(jj)+1,dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(numl(jj)+1,:),angleDropMP,Cang,WeSB,RvTent);
                if abs(errortan(3,jj+1))<abs(errortan(4,jj+1))
                    numlTent = numl(jj);
                    etaTent = etaprob(:,3);
                    phiTent = phiprob(:,3);
                    psNew = psprob(:,3);
                    zTent = zprob(3);
                    vzTent = vzprob(3);
                else
                    co = find(numl(jj:-1:1)~=numl(jj)+2,1);%I think I don't need this and I can just replace the first argument of solven by numl(jj)
                    [~,~,~,~,~,errortan(5,jj+1)] = ...
                        solvenDDCusp(numl(jj-co+1),numl(jj)+2,dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                        We,Ma,zs,IntMat(numl(jj)+2,:),angleDropMP,Cang,WeSB,RvTent);
                    if abs(errortan(4,jj+1)) < abs(errortan(5,jj+1))
                        numlTent = numl(jj)+1;
                        etaTent = etaprob(:,4);
                        phiTent = phiprob(:,4);
                        psNew = psprob(:,4);
                        zTent = zprob(4);
                        vzTent = vzprob(4);
                    else
                        tvec = [tvec(1:jj),tvec(jj)/2+tvec(jj+1)/2,tvec(jj+1:end)];
                        jj = jj-1; 
                        reduc = 1;
                    end
                end
            end
        elseif numl(jj) > nlmax(jj)-1.5 && numl(jj) < nlmaxTent-.5 %i.e. if last number of contacted points is nlmax-1
            co = find(numl(jj:-1:1)~=numl(jj),1);
            [etaprob(:,3),phiprob(:,3),zprob(3),vzprob(3),psprob(1:numl(jj),3),errortan(3,jj+1)] = ...    
                solvenDDCusp(numl(jj-co+1),numl(jj),dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                We,Ma,zs,IntMat(numl(jj),:),angleDropMP,Cang,WeSB,RvTent);
            co = find(numl(jj:-1:1)~=numl(jj)-1,1);
            [etaprob(:,2),phiprob(:,2),zprob(2),vzprob(2),psprob(1:numl(jj)-1,2),errortan(2,jj+1)] = ...    
                solvenDDCusp(numl(jj-co+1),numl(jj)-1,dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                We,Ma,zs,IntMat(numl(jj)-1,:),angleDropMP,Cang,WeSB,RvTent);
            if abs(errortan(2,jj+1))<abs(errortan(3,jj+1))
                co = find(numl(jj:-1:1)~=numl(jj)-2,1);
                [~,~,~,~,~,errortan(1,jj+1)] = ...    
                    solvenDDCusp(numl(jj-co+1),numl(jj)-2,dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(numl(jj)-2,:),angleDropMP,Cang,WeSB,RvTent);
                if abs(errortan(2,jj+1)) < abs(errortan(1,jj+1))
                    numlTent = numl(jj)-1;
                    etaTent = etaprob(:,2);
                    phiTent = phiprob(:,2);
                    psNew = psprob(:,2);
                    zTent = zprob(2);
                    vzTent = vzprob(2);
                else
                    tvec = [tvec(1:jj),tvec(jj)/2+tvec(jj+1)/2,tvec(jj+1:end)];
                    jj = jj-1; 
                    reduc = 1;
                end
            else
                co = find(numl(jj:-1:1)~=numl(jj)+1,1);
                [etaprob(:,4),phiprob(:,4),zprob(4),vzprob(4),psprob(1:numl(jj)+1,4),errortan(4,jj+1)] = ...    
                    solvenDDCusp(numl(jj-co+1),numl(jj)+1,dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(numl(jj)+1,:),angleDropMP,Cang,WeSB,RvTent);
                if abs(errortan(3,jj+1)) < abs(errortan(4,jj+1))
                    numlTent = numl(jj);
                    etaTent = etaprob(:,3);
                    phiTent = phiprob(:,3);
                    psNew = psprob(:,3);
                    zTent = zprob(3);
                    vzTent = vzprob(3);
                else
                    numlTent = numl(jj)+1;
                    etaTent = etaprob(:,4);
                    phiTent = phiprob(:,4);
                    psNew = psprob(:,4);
                    zTent = zprob(4);
                    vzTent = vzprob(4);
                end
            end
        elseif numl(jj) == nlmaxTent %i.e. if last number of contact points was nlmax
            co = find(numl(jj:-1:1)~=numl(jj),1);
            [etaprob(:,3),phiprob(:,3),zprob(3),vzprob(3),psprob(1:numl(jj),3),errortan(3,jj+1)] = ...    
                solvenDDCusp(numl(jj-co+1),numl(jj),dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                We,Ma,zs,IntMat(numl(jj),:),angleDropMP,Cang,WeSB,RvTent);
            co = find(numl(jj:-1:1)~=numl(jj)-1,1);
            [etaprob(:,2),phiprob(:,2),zprob(2),vzprob(2),psprob(1:numl(jj)-1,2),errortan(2,jj+1)] = ...    
                solvenDDCusp(numl(jj-co+1),numl(jj)-1,dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                We,Ma,zs,IntMat(numl(jj)-1,:),angleDropMP,Cang,WeSB,RvTent);
            if abs(errortan(2,jj+1)) < abs(errortan(3,jj+1))
                co = find(numl(jj:-1:1)~=numl(jj)-2,1);
                [~,~,~,~,~,errortan(1,jj+1)] = ...
                    solvenDDCusp(numl(jj-co+1),numl(jj)-2,dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(numl(jj)-2,:),angleDropMP,Cang,WeSB,RvTent);
                if abs(errortan(2,jj+1)) < abs(errortan(1,jj+1))
                    numlTent = numl(jj)-1;
                    etaTent = etaprob(:,2);
                    phiTent = phiprob(:,2);
                    psNew = psprob(:,2);
                    zTent = zprob(2);
                    vzTent = vzprob(2);
                else
                    tvec = [tvec(1:jj),tvec(jj)/2+tvec(jj+1)/2,tvec(jj+1:end)];
                    jj = jj-1; 
                    reduc = 1;
                end
            else
                numlTent = numl(jj);
                etaTent = etaprob(:,3);
                phiTent = phiprob(:,3);
                psNew = psprob(:,3);
                zTent = zprob(3);
                vzTent = vzprob(3);
            end
        else
            co = find(numl(jj:-1:1)~=numl(jj)-1,1);
            [etaprob(:,2),phiprob(:,2),zprob(2),vzprob(2),psprob(1:numl(jj)-1,2),errortan(2,jj+1)] = ...    
                solvenDDCusp(numl(jj-co+1),numl(jj)-1,dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                We,Ma,zs,IntMat(numl(jj)-1,:),angleDropMP,Cang,WeSB,RvTent);
            if abs(errortan(2,jj+1)) < 4
                co = find(numl(jj:-1:1)~=numl(jj)-2,1);
                [~,~,~,~,~,errortan(1,jj+1)] = ...
                    solvenDDCusp(numl(jj-co+1),numl(jj)-2,dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(numl(jj)-2,:),angleDropMP,Cang,WeSB,RvTent);
                if abs(errortan(2,jj+1)) < abs(errortan(1,jj+1))
                    numlTent = numl(jj)-1;
                    etaTent = etaprob(:,2);
                    phiTent = phiprob(:,2);
                    psNew = psprob(:,2);
                    zTent = zprob(2);
                    vzTent = vzprob(2);
                else
                    tvec = [tvec(1:jj),tvec(jj)/2+tvec(jj+1)/2,tvec(jj+1:end)];
                    jj = jj-1; 
                    reduc = 1;
                end
            else
                tvec = [tvec(1:jj),tvec(jj)/2+tvec(jj+1)/2,tvec(jj+1:end)];
                jj = jj-1; 
                reduc = 1;
            end
        end
        
        if ll == 100
            tvec = [tvec(1:jj),tvec(jj)/2+tvec(jj+1)/2,tvec(jj+1:end)];
            jj = jj-1;
            reduc = 1;
        end
        

        % verifying convergence of the pressure field
        if reduc == 0 %if there was no reduction of time step
            if norm(psNew,1) == 0
                %-%-B2New = 0;
                %-%-B3New = 0;
                B_l_ps_new = zeros(1, N);
            else
                nb_contact_points = nlmaxTent-find(flipud(psNew),1)+1;%Number of nodes in which the pressure needs to be integrated
                %against each harmonic 
                
                %%#---
                % We have that B_l = (2l+1)/2 * int(Ps * P_m * sin(theta))
                
                % the next line will project the pressur field onto th
                % spherical harmonics (i.e Legendre Polynomials)
%                 B_l_ps_new = arrayfun(@(m) find_harmonic_coefficient(thetaVec(1:(nb_contact_points+1)), ...
%                     psNew(1:nb_contact_points), m, nb_contact_points, LEGENDRE_POLYNOMIALS{m}), 1:N);
                
                f = @(thetas) interp1(thetaVec(1:(nb_contact_points+1)), [psNew(1:nb_contact_points)', 0], thetas, 'linear',  0); 
                endpoints = [thetaVec(nb_contact_points+1), thetaVec(1)];
                B_l_ps_new = project_amplitudes(f, N, endpoints, PROBLEM_CONSTANTS, true);   
            end
             
            
            %[amplitudes_new, amplitudes_velocities_new] = solve_ODE_unkown(N, dt, Ra, ...
            %    omegas_frequencies, B_l_ps_old, B_l_ps_new,  ODE_inverse_matrices, ...
            %    ODE_matrices, amplitudes_old, amplitudes_velocities_old);
            [amplitudes_new, velocities_new] = solve_ODE_unkown(nan, B_l_ps_new, dt, ...
                previous_conditions, PROBLEM_CONSTANTS);
            %#---
            
            nb_points = max(length(psTent),length(psNew));%number of points in which the pressure needs to be compared
            
            % #---
            err = norm([psTent;zeros(length(length(psTent)+1:nb_points),1);amplitudes_tent']-...
                       [psNew ;zeros(length(length(psNew )+1:nb_points),1);amplitudes_new'],1);
                   
            %-%-err = norm([psTent;zeros(length(length(psTent)+1:nb_points),1);A2Tent;A3Tent]-...
            %-%-           [psNew ;zeros(length(length(psNew )+1:nb_points),1);A2New ;A3New ],1);
            %-%-if norm([psNew;A2New;A3New],1) > 0
            %-%-    errorP = err/norm([psNew;A2New;A3New],1);
            if norm([psNew;amplitudes_new'],1) > 0
                errorP = err/norm([psNew;amplitudes_new'],1);    
            else
                errorP = err;
            end
            if errorP < tolP % Finally accept solution
                numl(jj+1) = numlTent;
                eta1 = etaTent;
                phi1 = phiTent;
                ps1 = psNew;
                z(jj+1) = zTent;
                vz(jj+1) = vzTent;
                %#---
                oscillation_amplitudes(:, jj + 1) = amplitudes_new;
                Rv(jj+1) = zs_from_spherical(pi, amplitudes_new);
                amplitudes_old = amplitudes_new;
                amplitudes_velocities_old = velocities_new;
                B_l_ps_old = B_l_ps_new;
                

                nlmax(jj+1) = nlmaxTent;
                etaOri(jj+1) = eta1(1);

                jj0 = floor(jj/nsteps);
                jj1 = round(jj-jj0*nsteps);
            %     if jj<10*nsteps || jj>=222*nsteps %this saves the start
            %     and the steady state in the periodic case
                etaMatPer(:,jj1+1) = eta1;
                phiMatPer(:,jj1+1) = phi1;
                psMatPer{jj1+1} = ps1;
                if jj1 == nsteps-1 
                    if runNumber == 0
                        tiempoComp(jj0+1)=toc(tstart);
                    end
                    save(['etaMatPer',num2str(jj0+1),'.mat'],'etaMatPer')
                    save(['phiMatPer',num2str(jj0+1),'.mat'],'phiMatPer')
                    save(['psMatPer',num2str(jj0+1),'.mat'],'psMatPer')

                    save('etaOri.mat','etaOri')
                    save('z.mat','z')
                    save('vz.mat','vz')
                    save('tvec.mat','tvec')
                    save('numl.mat','numl')
                    save('errortan.mat','errortan')
                    % s ave('oscillation_amplitudes.mat', 'oscillation_amplitudes');
                end

                if  zTent > 1.5 && numlTent == 0
                    tend = t;
                    tvec = tvec(1:(jj+1));
                    numl = numl(1:(jj+1));
                    z = z(1:(jj+1));
                    vz = vz(1:(jj+1));
                    oscillation_amplitudes = oscillation_amplitudes(:, 1:(jj + 1)); 
                    etaOri = etaOri(1:(jj+1));
                    etaMatPer = etaMatPer(:,1:(jj1+1));
                    phiMatPer = phiMatPer(:,1:(jj1+1));
                    Rv = Rv(1:(jj+1));
                end

                etao = eta1;
                phio = phi1;
                pso = ps1;
            
                xs = dr*(0:nlmax(jj+1)-1);
                zsplot = zs(1:nlmax(jj+1))+RvTent+z(jj+1);
                plot([-fliplr(xs(2:end)),xs],[flipud(zsplot(2:end));zsplot],'k','Linewidth',2);
                hold on
                thetaplot = linspace(0, thetaVec(end), 200);%-%-0:thetaVec(end)/400:thetaVec(end);
                %-%-xsTop = xsoftheta(thetaplot,A2New,A3New);
                %-%-zsTop = zsoftheta(thetaplot,A2New,A3New);
                zsTop = zs_from_spherical(thetaplot, amplitudes_new);
                xsTop = r_from_spherical(thetaplot, amplitudes_new); 
                plot([-xsTop(end:-1:2), xsTop],[zsTop(end:-1:2), zsTop]+zTent,'k','Linewidth',2);
                width = round(nr * 6/D);
                plot([-fliplr(xplot(2:width)),xplot(1:width)],[flipud(eta1(2:width));eta1(1:width)],'LineWidth',2);
                hold off
                axis equal
                title(['   t = ',num2str(t),'   ','nl = ',num2str(numl(jj+1))],'FontSize',16);
                grid on
                %set(gca,'xlim',[-6 6])
                drawnow;
            else
                
                %i.e. didn't attain convergence yet
                psTent = psNew;
                
                %#---
                amplitudes_tent = amplitudes_new;
                B_l_ps_tent = B_l_ps_new;
                %#---

                RmaxTent = r_from_spherical(maximum_contact_radius(amplitudes_tent), amplitudes_tent);
                
                %RmaxTent = r_from_spherical(theta_max_radius(amplitudes_tent), amplitudes_tent);
                %-%-RmaxTent = xsoftheta(mod(abs(thetaMax(A2Tent,A3Tent)),pi),A2Tent,A3Tent);
                nlmaxTent = floor(RmaxTent/dr)+1;
                
                thetaVec = theta_from_cylindrical(dr*(0:(nlmaxTent-1)), oscillation_amplitudes(:, jj));
                 
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

%                 thetaVecMP = zeros(nlmaxTent,1);
%                 thetaVecMP(1) = mod(abs(thetaofxs(.5*dr,A2Tent,A3Tent,pi)),pi);
%                 for ii = 2:nlmaxTent-1
%                     thetaVecMP(ii) = mod(abs(thetaofxs(dr*(ii-.5),A2Tent,A3Tent,thetaVecMP(ii-1))),pi);
%                 end
%                 thetaVecMP(nlmaxTent) = mod(abs(thetaofxs(dr*(nlmaxTent-1),A2Tent,A3Tent,thetaVecMP(nlmaxTent-1))),pi);
%                 ysp = ysprimeoftheta(thetaVecMP,A2Tent,A3Tent);
%                 xsp = xsprimeoftheta(thetaVecMP,A2Tent,A3Tent);
%                 angleDropMP = zeros(nlmaxTent,1);
%                 if xsp(end) == 0
%                     tanDrop(1:nlmaxTent-1) = ysp(1:end-1)./xsp(1:end-1);
%                     angleDropMP(1:nlmaxTent-1) = atan(tanDrop(1:nlmaxTent-1));
%                     angleDropMP(nlmaxTent) = pi/2;
%                 else
%                     tanDrop(1:nlmaxTent) = ysp./xsp;
%                     angleDropMP(1:nlmaxTent) = atan(tanDrop(1:nlmaxTent));
%                 end
                psprob = zeros(nlmaxTent,5);%zeroing the vector of potential pressures
            end
        else
            if 1/(dt * nsteps) >= 2^20
                warning("Step size has been made too small (%.3e). Stopped the execution of the program", dt);
                t = inf;
            end
        end
    end
    

    tComp = toc(tstart);
    if jj1==0 && tComp > tmax %#---
        runNumber = runNumber+1;
        save('runNumber.mat','runNumber')
        tstop = t;
        save(['tstop',num2str(runNumber),'.mat'],'tstop')
        jjstop = jj;
        save(['jjstop',num2str(runNumber),'.mat'],'jjstop')
        save(['etao',num2str(runNumber),'.mat'],'etao')
        save(['phio',num2str(runNumber),'.mat'],'phio')
        save(['pso',num2str(runNumber),'.mat'],'pso')

        zrestart = z(jj+1);
        vzrestart = vz(jj+1);
        trestart = tvec(jj+1);
        numlrestart = numl(jj+1);
        save(['zrestart',num2str(runNumber),'.mat'],'zrestart')
        save(['vzrestart',num2str(runNumber),'.mat'],'vzrestart')
        save(['trestart',num2str(runNumber),'.mat'],'trestart')
        save(['numlrestart',num2str(runNumber),'.mat'],'numlrestart')

        save('etaOri.mat','etaOri')
        %-%-save('A2.mat','A2')
        %-%-save('V2.mat','V2')
        %-%-save('A3.mat','A3')
        %-%-save('V3.mat','V3')
        save('z.mat','z')
        save('vz.mat','vz')
        save('tvec.mat','tvec')
        save('numl.mat','numl')
        save('nlmax.mat','nlmax')
        save('errortan.mat','errortan')
        if runNumber == 1
            save('tiempoComp.mat','tiempoComp')
        end
        
        %exit
    end
end

runNumber = runNumber+1;
tstop = t;
save(['tstop',num2str(runNumber),'.mat'],'tstop')
jjstop = jj;
save(['jjstop',num2str(runNumber),'.mat'],'jjstop')
save(['etao',num2str(runNumber),'.mat'],'etao')
save(['phio',num2str(runNumber),'.mat'],'phio')
save(['pso',num2str(runNumber),'.mat'],'pso')

zrestart = z(jj+1);
vzrestart = vz(jj+1);
trestart = tvec(jj+1);
numlrestart = numl(jj+1);
save(['zrestart',num2str(runNumber),'.mat'],'zrestart')
save(['vzrestart',num2str(runNumber),'.mat'],'vzrestart')
save(['trestart',num2str(runNumber),'.mat'],'trestart')
save(['numlrestart',num2str(runNumber),'.mat'],'numlrestart')

save('etaOri.mat','etaOri')
save('z.mat','z')

save('vz.mat','vz')
save('tvec.mat','tvec');
save('nlmax.mat','nlmax');
save('numl.mat','numl');
save('errortan.mat','errortan');
save('oscillation_amplitudes.mat', 'oscillation_amplitudes');
save('Rv.mat', 'Rv');
runNumber = -10;
save('runNumber.mat','runNumber');


% 
% function [amplitudes_tent, amplitudes_velocities_tent] = solve_EDO(N, dt, Ra,...
%     omegas_frequencies, B_l_ps_old, B_l_ps_tent, ODE_inverse_matrices, ODE_matrices, amplitudes_old, amplitudes_velocities_old)
% %     if norm(B_l_ps_old) + norm(B_l_ps_tent) == 0
% %         
% %        return 
% %     end
%     
%     Y_old = zeros(2, N);
%     C_tent= zeros(2, N);
%     C_old = zeros(2, N);
%     Y_tent= zeros(2, N);
%     amplitudes_tent = zeros(1, N);
%     amplitudes_velocities_tent = zeros(1, N);
%     
%     for ii = 2:N % We dont care about ii = 1!
%         Y_old(:, ii) = ODE_inverse_matrices(:, :, ii) * [amplitudes_old(ii); amplitudes_velocities_old(ii)];
%         C_old(:, ii) = ODE_inverse_matrices(:, :, ii) * [0; -ii*B_l_ps_old(ii)/Ra];
%         C_tent(:, ii) = ODE_inverse_matrices(:, :, ii)* [0; -ii*B_l_ps_tent(ii)/Ra];
%         
%         %Y_tent(:, ii) = (eye(2)-dt*diag([1.0i * omegas_frequencies(ii), -1.0i * omegas_frequencies(ii)]))\(dt*C_tent(:, ii) + Y_old(:, ii));
%         Y_tent(:, ii) = diag(exp([1.0i * dt * omegas_frequencies(ii), -1.0i * dt* omegas_frequencies(ii)])) * ...
%             (Y_old(:, ii) + dt/2 * C_old(:, ii)) + dt/2 * C_tent(:, ii);
%         X = ODE_matrices(:, :, ii) * Y_tent(:, ii);
%         amplitudes_tent(ii) = X(1);
%         amplitudes_velocities_tent(ii) = X(2);
%     end
% end
% 
% function I = find_harmonic_coefficient(angles, pressure_vector, m, nb_contact_points, legendrePol)
%     % This function will return the coefficient corresponding to the
%     % expansion in Legendre POlynomials of the pressure vector,
%     % corresponding to the mth Legendre Polynomial.
%     
%     if m < 2
%         I = 0;
%     else
%         if size(pressure_vector, 1) > 1; pressure_vector = pressure_vector'; end
%         pressure_extended = [pressure_vector(1:nb_contact_points), -pressure_vector(nb_contact_points)];
%         y = @(angle) interp1(angles, pressure_extended, angle) .* ...
%                 legendrePol(cos(angle)) .* sin(angle); % TODO: Interpolate assuming linearity on radius
%         I = (m + 1/2) * integral(y, (angles(nb_contact_points) + angles(nb_contact_points+1))/2, pi, "RelTol", 1e-4);
%     end
% end