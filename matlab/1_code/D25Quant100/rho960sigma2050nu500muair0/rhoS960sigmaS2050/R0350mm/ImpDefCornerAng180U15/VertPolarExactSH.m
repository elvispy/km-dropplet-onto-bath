%clear
close all
%clc

tstart = tic;
tic
%data in cgs
tmax = 100;


if exist('z.mat', 'file') == 2
   % error("Exporting data is going to be overwritten. Please re-allocate files to avoid loss of data");
end


U0 = 15; %impact velocity in cm/s (unit of velocity for the problem)
Ang = 180; %contact angle to be imposed

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
N = 15; % Number of harmonics contributing to the oscillation
% #---0

%Unit of time
% T = Ro/U0; %base time is seconds
X = Ro; % Characteristic length
T = sqrt(rhoS * Ro^3/sigmaS); % Characteristic time

%Dimensionless numbers that depend on U0
Dr = rhoS/rho; Sr = sigmaS/sigma;
Re = Ro^2/(nu * T); % Ro*U0/nu; 
Fr = sigma/(g * rho * Ro^2);% U0^2/(g*Ro); 
We = Sr * rhoS * Ro.^3 / (sigmaS * T^2); % rho*Ro*U0^2/sigma; 
%Dr = We; %rhoS*Ro*U0^2/sigma;
WeS  = We;%rhoS*Ro*U0^2/sigmaS; %This name may not be the best, the surface tension is that of the 
%bath at least in one place
Cang = (Ang/180)*pi; %contact angle to be imposed

%Physical parameter
tend = 9; %Earliest possible end of simulation in characteristic units

%Inintial conditions for the fluid
t = 0;
etao = zeros(nr,1); %initial surface elevation
phio = zeros(nr,1); %initial surface potential

%Numerical Simulation parameters
nsteps = 100; %minimum number of timesteps in one unit of time
dtb = 1/nsteps; %basic timestep (gets halved as needed over impacts)
steps = ceil((tend-t)/dtb); %estimated minimum number of timesteps

%Zeroing result storing variables
etaOri = zeros(1,steps+1);%height of the surface below the south pole
z = zeros(1,steps+1);%height of the centre of mass
vz = zeros(1,steps+1);%speed of the centre of mass
numl = zeros(1,steps+1);%number of pressed mesh points at each time step
tvec = t:(dtb):tend+1; tvecOri = tvec;%vector of times assuming no refinement has happened
%plus some extra time just in case the simulation needs to run longer
% #--- 
dt = tvec(2) - tvec(1); indexes_to_save = zeros(steps + 1, 1);
current_to_save = 2; indexes_to_save(1) = 1;
oscillation_amplitudes = zeros(N, steps + 1); % Variable to store
Rv = -ones(1, steps+1);
% the time dependent amplitude of all the SH
oscillation_velocities = zeros(N, steps+1);
% #--- 
nlmax = zeros(1,steps+1);%Variable to store the number of nodes spanned by the deformed droplet

tolP = 1E-6; %error tolerance for the pressure field and deformation 

save('ProblemConditions.mat', "T", "N", "U0", "Ang", "Re", "Fr", "We", ...
"WeS", "Cang", "tend", "nsteps", "dtb" );

%Drop oscillation frequencies
% #--- 
f = @(n) sqrt(n.*(n+2).*(n-1)./1); % Changed frequency
omegas_frequencies = f(1:N)';

% #---oscillation_amplitudes = zeros(N, steps + 1);
amplitudes_old = oscillation_amplitudes(:, 1);
amplitudes_velocities_old = oscillation_velocities(:, 1);
B_l_ps_old = zeros(1, N);

% # ---
z(1) = -1* zs_from_spherical(pi, oscillation_amplitudes(:, 1));% -1*zsoftheta(pi,A2(1),A3(1)); %height of the centre of mass (CoM) in dimensionless units,

% zsoftheta(pi,A2(1),A3(1)) gives the height of the south pole with
% respect to the CoM, z(1) is chosen so that the drop is just about to touch down
vz(1) = -abs(U0/ (X/T)); %Initial velocity of the CoM in dimesionless units


current_conditions = struct("deformation_amplitudes", amplitudes_old, ...
    "deformation_velocities", amplitudes_velocities_old, ...
    "pressure_amplitudes", B_l_ps_old, "dt", dt, "nb_harmonics", N,  ...
    "current_time", 0, ...
    "center_of_mass", z(1), "center_of_mass_velocity", vz(1), ...
    "nb_contact_points", 0);

previous_conditions = {current_conditions, current_conditions}; 

% f = @(n)  sqrt(n .* (n+2) .* (n-1) / WeS);
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

tentative_index = 0;%iteration counter

errortan = zeros(5,steps+1);%tangency error recorder

ps_accepted = [];


%If there were some initial pressure acting on the surface and sphere I
%would have to change this bit here to reflect the presure distribution

%zeroing variable that records each part of the sequence of surface states
etaMatPer = zeros(length(etao),nsteps); 
etas      = zeros(length(etao), steps); 

phiMatPer = zeros(length(phio),nsteps);
psMatPer = cell(1,nsteps); % ??? why cell
%Storing initial surface state
etaMatPer(:,1) = etao;
phiMatPer(:,1) = phio;
psMatPer{1} = zeros(quant+1,1);

%zeroing the ceiling functions
zs = zeros(nr,1);

jj1 = 1; %partial results savings  counter

PROBLEM_CONSTANTS = struct("froude_nb", Fr, "weber_nb", We, ...
    "nb_harmonics", N, ...
    "density_ratio", Dr, ...
    "omegas_frequencies", omegas_frequencies, ...
    "spatial_tol", dr, ...
    "DEBUG_FLAG", true, "linear_on_theta", true, ...
    "Ra", Ra, "interpolation_number", 10);
                            %"pressure_unit", pressure_unit, ...
                            %"CM", 9, ...
                            %"PG", 2, ...
                            %"KILL_OUTSIDE", true, ...
                            %"wigner3j", {precomputed_wigner(harmonics_qtt)}, ...
        

%% Main Loop
while (t<tend) %#-- || jj1>.5) 
%     if toc > 120
%         break;
%     end
    tentative_index = tentative_index+1;
    t = tvec(tentative_index+1);
    dt = t - tvec(tentative_index);
    %zeroing the tentative solution variables
    etaprob = zeros(nr,5);
    phiprob = zeros(nr,5);
    vzprob  = zeros(1,5);
    zprob   = zeros(1,5);
    errortan(:,tentative_index+1) = 4*ones(5,1);
    
    psTent = ps_accepted; %Tentative pressure distribution (we start with the previous pressure)
    
    RmaxOld = r_from_spherical(maximum_contact_radius(oscillation_amplitudes(:, tentative_index)), oscillation_amplitudes(:, tentative_index));
  
    %(i.e. where the tangent plane to the droplet is vertical)
    nlmax(tentative_index) = floor(RmaxOld/dr)+1;%max number of contact points
    
    thetaVec = theta_from_cylindrical(dr*(0:(nlmax(tentative_index)-1)), oscillation_amplitudes(:, tentative_index)); % zeros(1,nlmax(jj));%initialising vector of angles of pressed positions

    % Spherical Harmonics modes
    if norm(psTent,1) == 0
        %B2Tent = 0;
        %B3Tent = 0;
        B_l_ps_tent = zeros(1, N);
    else
        nb_contact_points = nlmax(tentative_index)-find(flipud(psTent),1)+1; %Number of nodes contact points%
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
                angles = [thetaVec(1:(nb_contact_points)), (2 * thetaVec(nb_contact_points) - thetaVec(nb_contact_points-1))];
            else
                angles = thetaVec(1:(nb_contact_points+1));
            end
            f = @(thetas) interp1(angles, [psTent(1:nb_contact_points)', 0], thetas, 'linear',  0); 
            endpoints = [angles(end), angles(1)];
            B_l_ps_tent = project_amplitudes(f, N, endpoints, PROBLEM_CONSTANTS, true); 
        end
    end    

    [amplitudes_tent, velocities_tent] = solve_ODE_unkown(nan, B_l_ps_tent, dt, ...
        previous_conditions, PROBLEM_CONSTANTS);

    
    RmaxTent = r_from_spherical(maximum_contact_radius(oscillation_amplitudes(:, tentative_index)), oscillation_amplitudes(:, tentative_index));
  
    nlmaxTent = floor(RmaxTent/dr)+1;
    thetaVec  = theta_from_cylindrical(dr*(0:(nlmaxTent-1)), oscillation_amplitudes(:, tentative_index)); % zeros(1,nlmaxTent);

    RvTent = zs_from_spherical(pi, amplitudes_tent);
    zs(1:nlmaxTent) = zs_from_spherical(thetaVec, amplitudes_tent)' - RvTent; %TODO: Check that matrix dimensions agree.
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
                solveDD0(dt,z(tentative_index),vz(tentative_index),etao,phio,nr,Re,Delta,DTN,Fr,We,zs,RvTent);
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
                    solvenDDCusp(numl(tentative_index-co+1),1,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(1,:),angleDropMP,Cang,Dr,RvTent);
                co = find(numl(tentative_index:-1:1)~=2,1);
                [~,~,~,~,~,errortan(5,tentative_index+1)] = ...    
                    solvenDDCusp(numl(tentative_index-co+1),2,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
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
                solveDD0(dt,z(tentative_index),vz(tentative_index),etao,phio,nr,Re,Delta,DTN,Fr,We,zs,RvTent);
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
                    solvenDDCusp(numl(tentative_index-co+1),1,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(1,:),angleDropMP,Cang,Dr,RvTent);
                co = find(numl(tentative_index:-1:1)~=2,1);
                [etaprob(:,4),phiprob(:,4),zprob(4),vzprob(4),psprob(1:2,4),errortan(4,tentative_index+1)] = ...    
                    solvenDDCusp(numl(tentative_index-co+1),2,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
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
                        solvenDDCusp(numl(tentative_index-co+1),3,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
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
                solveDD0(dt,z(tentative_index),vz(tentative_index),etao,phio,nr,Re,Delta,DTN,Fr,We,zs,RvTent);
            if abs(errortan(1,tentative_index+1))<.5
                tvec = [tvec(1:tentative_index),tvec(tentative_index)/2+tvec(tentative_index+1)/2,tvec(tentative_index+1:end)];
                tentative_index = tentative_index-1; 
                reduc = 1;
            else
                co = find(numl(tentative_index:-1:1)~=2,1);
                [etaprob(:,3),phiprob(:,3),zprob(3),vzprob(3),psprob(1:2,3),errortan(3,tentative_index+1)] = ...    
                    solvenDDCusp(numl(tentative_index-co+1),2,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(2,:),angleDropMP,Cang,Dr,RvTent);
                co = find(numl(tentative_index:-1:1)~=1,1);
                [etaprob(:,2),phiprob(:,2),zprob(2),vzprob(2),psprob(1,2),errortan(2,tentative_index+1)] = ...    
                    solvenDDCusp(numl(tentative_index-co+1),1,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
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
                        solvenDDCusp(numl(tentative_index-co+1),3,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
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
                            solvenDDCusp(numl(tentative_index-co+1),4,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
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
                solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index),dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                We,Ma,zs,IntMat(numl(tentative_index),:),angleDropMP,Cang,Dr,RvTent);
            co = find(numl(tentative_index:-1:1)~=numl(tentative_index)-1,1);
            [etaprob(:,2),phiprob(:,2),zprob(2),vzprob(2),psprob(1:numl(tentative_index)-1,2),errortan(2,tentative_index+1)] = ...    
                solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index)-1,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                We,Ma,zs,IntMat(numl(tentative_index)-1,:),angleDropMP,Cang,Dr,RvTent);
            if abs(errortan(2,tentative_index+1)) < abs(errortan(3,tentative_index+1))
                co = find(numl(tentative_index:-1:1)~=numl(tentative_index)-2,1);
                [~,~,~,~,~,errortan(1,tentative_index+1)] = ...
                    solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index)-2,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
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
                    solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index)+1,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
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
                        solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index)+2,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
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
                solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index),dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                We,Ma,zs,IntMat(numl(tentative_index),:),angleDropMP,Cang,Dr,RvTent);
            co = find(numl(tentative_index:-1:1)~=numl(tentative_index)-1,1);
            [etaprob(:,2),phiprob(:,2),zprob(2),vzprob(2),psprob(1:numl(tentative_index)-1,2),errortan(2,tentative_index+1)] = ...    
                solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index)-1,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                We,Ma,zs,IntMat(numl(tentative_index)-1,:),angleDropMP,Cang,Dr,RvTent);
            if abs(errortan(2,tentative_index+1))<abs(errortan(3,tentative_index+1))
                co = find(numl(tentative_index:-1:1)~=numl(tentative_index)-2,1);
                [~,~,~,~,~,errortan(1,tentative_index+1)] = ...    
                    solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index)-2,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
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
                    solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index)+1,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
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
                solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index),dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                We,Ma,zs,IntMat(numl(tentative_index),:),angleDropMP,Cang,Dr,RvTent);
            co = find(numl(tentative_index:-1:1)~=numl(tentative_index)-1,1);
            [etaprob(:,2),phiprob(:,2),zprob(2),vzprob(2),psprob(1:numl(tentative_index)-1,2),errortan(2,tentative_index+1)] = ...    
                solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index)-1,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                We,Ma,zs,IntMat(numl(tentative_index)-1,:),angleDropMP,Cang,Dr,RvTent);
            if abs(errortan(2,tentative_index+1)) < abs(errortan(3,tentative_index+1))
                co = find(numl(tentative_index:-1:1)~=numl(tentative_index)-2,1);
                [~,~,~,~,~,errortan(1,tentative_index+1)] = ...
                    solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index)-2,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
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
                solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index)-1,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                We,Ma,zs,IntMat(numl(tentative_index)-1,:),angleDropMP,Cang,Dr,RvTent);
            if abs(errortan(2,tentative_index+1)) < 4
                co = find(numl(tentative_index:-1:1)~=numl(tentative_index)-2,1);
                [~,~,~,~,~,errortan(1,tentative_index+1)] = ...
                    solvenDDCusp(numl(tentative_index-co+1),numl(tentative_index)-2,dt,z(tentative_index),vz(tentative_index),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
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
        
        if ll == 100
            tvec = [tvec(1:tentative_index),tvec(tentative_index)/2+tvec(tentative_index+1)/2,tvec(tentative_index+1:end)];
            tentative_index = tentative_index-1;
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
            
            nb_points = max(length(psTent),length(psNew));%number of points in which the pressure needs to be compared
            
            err = norm([psTent;zeros(length(length(psTent)+1:nb_points),1);amplitudes_tent']-...
                       [psNew ;zeros(length(length(psNew )+1:nb_points),1);amplitudes_new'],1);
                   
            if norm([psNew;amplitudes_new'],1) > 0
                errorP = err/norm([psNew;amplitudes_new'],1);    
            else
                errorP = err;
            end
            if errorP < tolP % Finally accept solution
                numl(tentative_index+1) = numlTent;
                eta_accepted = etaTent;
                phi_accepted = phiTent;
                ps_accepted = psNew;
                z(tentative_index+1) = zTent;
                vz(tentative_index+1) = vzTent;
                %#---
                oscillation_amplitudes(:, tentative_index + 1) = amplitudes_new;
                Rv(tentative_index+1) = zs_from_spherical(pi, amplitudes_new);
                amplitudes_old = amplitudes_new;
                amplitudes_velocities_old = velocities_new;
                B_l_ps_old = B_l_ps_new;
                
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
                %jj0 = floor(tentative_index/nsteps);
                %jj1 = round(tentative_index-jj0*nsteps);
            %     if jj<10*nsteps || jj>=222*nsteps %this saves the start
            %     and the steady state in the periodic case
                if t >= tvecOri(current_to_save)
                    etaMatPer(:,current_to_save) = eta_accepted;
                    phiMatPer(:,current_to_save) = phi_accepted;
                    psMatPer{current_to_save} = ps_accepted;
                    indexes_to_save(current_to_save) = tentative_index;
                    current_to_save = current_to_save + 1;
                end
%                 if jj1 == nsteps-1 
%                     if runNumber == 0
%                         tiempoComp(jj0+1)=toc(tstart);
%                     end
%                     save(['etaMatPer',num2str(jj0+1),'.mat'],'etaMatPer')
%                     save(['phiMatPer',num2str(jj0+1),'.mat'],'phiMatPer')
%                     save(['psMatPer',num2str(jj0+1),'.mat'],'psMatPer')
% 
%                     save('etaOri.mat','etaOri')
%                     save('z.mat','z')
%                     save('vz.mat','vz')
%                     save('tvec.mat','tvec')
%                     save('numl.mat','numl')
%                     save('errortan.mat','errortan')
%                     % s ave('oscillation_amplitudes.mat', 'oscillation_amplitudes');
%                 end

                if zTent > 1.5 && numlTent == 0
                    tend = t;
                end

                etao = eta_accepted;
                phio = phi_accepted;
                pso = ps_accepted;
            
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
                title(['   t = ',num2str(t),'   ','nl = ',num2str(numl(tentative_index+1))],'FontSize',16);
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
        else
            if dt < 1e-15
                t = inf;
                error("Step size has been made too small (%.3e). Stopped the execution of the program", dt);
                
            end
        end %it time step was not reduced
    end % inner while   
end % Outer while

% runNumber = runNumber+1;
tstop = t;

if t < inf
    indexes_to_save = indexes_to_save(1:(current_to_save-1));
    z = z(indexes_to_save); save('z.mat','z')
    etaOri = etaOri(indexes_to_save); save('etaOri.mat','etaOri')
    etas = etas(:, indexes_to_save); save('etas.mat', 'etas');
    etaMatPer = etaMatPer(:,  1:(current_to_save-1)); save('etaMatPer.mat', 'etaMatPer');
    phiMatPer = phiMatPer(:,  1:(current_to_save-1)); save('phiMatPer.mat','phiMatPer');
    psMatPer  = psMatPer{1:(current_to_save-1)};  save('psMatPer.mat','psMatPer');

    vz = vz(indexes_to_save); save('vz.mat','vz');
    tvecOri = tvecOri(1:(current_to_save-1)); tvec = tvecOri; save('tvec.mat','tvec'); 
    nlmax = nlmax(indexes_to_save); save('nlmax.mat','nlmax');
    numl = numl(indexes_to_save); save('numl.mat','numl');
    % errrortan = errortan(indexes_to_save); save('errortan.mat','errortan');
    oscillation_amplitudes = oscillation_amplitudes(:, indexes_to_save); save('oscillation_amplitudes.mat', 'oscillation_amplitudes');
    Rv = Rv(indexes_to_save); save('Rv.mat', 'Rv');

    % 
    % save('etaOri.mat','etaOri')
    % save('z.mat','z')
    % 
    % save('vz.mat','vz')
    % save('tvec.mat','tvec');
    % save('nlmax.mat','nlmax');
    % save('numl.mat','numl');
    % save('errortan.mat','errortan');
    % save('oscillation_amplitudes.mat', 'oscillation_amplitudes');
    % save('Rv.mat', 'Rv');
    % 
end
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