%clear
close all
%clc

%Reset warning
lastwarn('', '');

tstart = tic;
tic
%data in cgs
tmax = 100;


if exist('z.mat', 'file') == 2
   % error("Exporting data is going to be overwritten.
end

if exist('oscillation_amplitudes.mat', 'file') == 2
   % error("Exporting data is going to be overwritten. Please re-allocate files to avoid loss of data");
end


U0 = 17; %impact velocity in cm/s (unit of velocity for the problem)
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
xplot = dr*(0:nr-1); 

cd(['rho',num2str(1000*rho),'sigma',num2str(round(100*sigma)),'nu',num2str(round(10000*nu)),'muair',num2str(muair)])

cd(['RhoS',num2str(rhoS*1000),'SigmaS',num2str(round(100*sigmaS))])
load('Ma.mat','Ma')%Dimensionless mass of sphere
load('Ra.mat','Ra')%Density ratio

cd(sprintf('R0%gmm',Ro*10000))

cd(['ImpDefCornerAng',num2str(Ang),'U',num2str(U0)])

tiempoComp = zeros(1,10); %just to check how long it takes to solve the first ten saving intervals

% #--- 
N = 10; % Number of harmonics contributing to the oscillation
% #---0

%Characteristic Units
L_unit = Ro; 
M_unit = rhoS * L_unit^3;
T = sqrt(rhoS * Ro^3/sigmaS); % Characteristic time
%    T = Ro/U0; 
T_unit = T;
V_unit = L_unit/T_unit;

%Dimensionless numbers for equations
Dr = rhoS/rho; %Sr = sigmaS/sigma;
Re = L_unit^2/(nu*T_unit);
Fr = L_unit/(g * T_unit^2);
We = rho * L_unit.^3 / (sigma * T_unit^2); 
WeS  = rhoS*Ro^3/(sigmaS * T_unit^2); %This is for the bath/dropplet interaction.

Cang = (Ang/180)*pi; %contact angle to be imposed

%Physical parameters
tend = 8; %Earliest possible end of simulation in characteristic units

%Inintial conditions for the fluid
t = 0;
etao = zeros(nr,1); %initial surface elevation
phio = zeros(nr,1); %initial surface potential

%Numerical Simulation parameters
nsteps = ceil(20/(2*pi) * N^(3/2)); %minimum number of timesteps in one unit of time
%nsteps = 100; %nsteps + 100 - mod(nsteps, 100);
dtb = 1/nsteps; %basic timestep (gets halved as needed over impacts)
steps = ceil((tend-t)/dtb); %estimated minimum number of timesteps

%Zeroing result storing variables
etaOri = zeros(1,steps+1);%height of the surface below the south pole
z = zeros(1,steps+1);%height of the centre of mass
vz = zeros(1,steps+1);%speed of the centre of mass
numl = zeros(1,steps+1);%number of pressed mesh points at each time step
tvec = zeros(1, steps+1); %vector of times assuming no refinement has happened
%plus some extra time just in case the simulation needs to run longer

dt = dtb; time_index = 0; save_index = 0; resetter = 0;
indexes_to_save = zeros(steps + 1, 1); indexes_to_save(1) = 1;
current_to_save = 2;
oscillation_amplitudes = zeros(N, steps + 1); % Harmonic amplitudes
Rv = -ones(1, steps+1);
% the time dependent amplitude of all the SH
oscillation_velocities = zeros(N, steps+1);
nlmax = zeros(1,steps+1);%Variable to store the number of nodes spanned by the deformed droplet

tolP = 1E-4; %error tolerance for the pressure field and deformation 
    
save('ProblemConditions.mat', "T", "N", "U0", "Ang", "Re", "Fr", "We", ...
    "WeS", "Cang", "tend", "nsteps", "dtb", "L_unit", "T_unit");

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

%zeroing variable that records each part of the sequence of surface states
etaMatPer = zeros(length(etao),nsteps); 
etas      = zeros(length(etao), steps);
phiMatPer = zeros(length(phio),nsteps);
psMatPer = cell(1,nsteps); % ??? why cell
%Storing initial surface state
etas(:, 1) = etao;
etaMatPer(:,1) = etao;
phiMatPer(:,1) = phio;
psMatPer{1} = zeros(quant+1,1);

%zeroing the ceiling functions
zs = zeros(nr,1);

jj1 = 1; %partial results savings  counter

PROBLEM_CONSTANTS = struct("We", We, "Ma", Ma, ...
    "Cang", Cang, "Dr", Dr, "Delta", Delta, ...
    "IntMat", IntMat, "nr", nr, "dr", dr, "DTN", DTN, "Fr", Fr, ...
    "froude_nb", Fr, "weber_nb", WeS, ...
    "nb_harmonics", N, "initial_dt", dtb, ...
    "omegas_frequencies", omegas_frequencies, ...
    "spatial_tol", dr, ...
    "DEBUG_FLAG", true, "linear_on_theta", true, ...
    "interpolation_number", 10, ...
    "Ra", Ra,  "Re", Re);
exit = false;
%% Main Loop
while tentative_index <= steps && t < tend && exit == false 
%     if toc > 120
%         break;
%     end
    if PROBLEM_CONSTANTS.DEBUG_FLAG == true
       fprintf("Outside %0.4g, %0.3e\n", t, dt); 
    end
    tentative_index = tentative_index+1;
    
    
    %zeroing the tentative solution variables
    etaprob = zeros(nr,5);
    phiprob = zeros(nr,5);
    vzprob  = zeros(1,5);
    zprob   = zeros(1,5);
    errortan(:,tentative_index+1) = Inf*ones(5,1);
    
    psTent = ps_accepted; %Tentative pressure distribution (we start with the previous pressure)
    
    %Shape at the start of the time step
    %x that corresponds to the max pressed radius 
    
    RmaxOld = r_from_spherical(maximum_contact_radius(oscillation_amplitudes(:, tentative_index)), oscillation_amplitudes(:, tentative_index));
  
    %(i.e. where the tangent plane to the droplet is vertical)
    nlmax(tentative_index) = floor(RmaxOld/dr)+1;%max number of contact points

    thetaVec = theta_from_cylindrical(dr*(0:(nlmax(tentative_index)-1)), oscillation_amplitudes(:, tentative_index)); % zeros(1,nlmax(jj));%initialising vector of angles of pressed positions

    %Spherical Harmonic modes
    if norm(psTent,1) == 0

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
    
    
    %#--- Solving the ODE

    [amplitudes_tent, velocities_tent] = solve_ODE_unkown(nan, B_l_ps_tent, dt, ...
        previous_conditions, PROBLEM_CONSTANTS);
    
    RmaxTent = r_from_spherical(maximum_contact_radius(oscillation_amplitudes(:, tentative_index)), oscillation_amplitudes(:, tentative_index));

    nlmaxTent = floor(RmaxTent/dr)+1;
    thetaVec  = theta_from_cylindrical(dr*(0:(nlmaxTent-1)), oscillation_amplitudes(:, tentative_index)); % zeros(1,nlmaxTent);

    RvTent = zs_from_spherical(pi, amplitudes_tent);
    zs(1:nlmaxTent) = zs_from_spherical(thetaVec, amplitudes_tent)' - RvTent;
    zs((nlmaxTent+1):nr) = Inf; 
    
    tanDrop = calculate_tan( dr * (1:nlmaxTent) - dr/2, amplitudes_tent)';
    angleDropMP(1:(nlmaxTent)) = atan(tanDrop(1:(nlmaxTent)));
    
    %% Inner while: Solving rigid problem
    psprob = zeros(nlmaxTent,5);%zeroing the vector of potential pressures
    errorP = 1; %error in the pressure field and amplitude of modes
    reduc = 0; %indicator of whether there was a reduction in the time-step size or not
    ll = 0; % Limiting while loop to 100 iterations
    while abs(errorP)>=tolP && reduc == 0 
        ll = ll + 1;
        
        if PROBLEM_CONSTANTS.DEBUG_FLAG == true
           %fprintf("Inside ll: %0.2g, psTent: %0.5g\n", ll, norm(psTent));
           %disp(errortan(:, tentative_index+1)');
        end
        

        [etaprob(:,3),phiprob(:,3),zprob(3),vzprob(3),psprob(1:numl(tentative_index),3),errortan(3,tentative_index+1)] = ...
            getNextStep(numl(tentative_index),numl(tentative_index),dt,z(tentative_index),vz(tentative_index),etaMatPer(:, tentative_index), phiMatPer(:, tentative_index), ...
            zs,RvTent, angleDropMP, PROBLEM_CONSTANTS, nlmaxTent);        
        if abs(errortan(3, tentative_index+1)) < 1e-5
            numlTent = numl(tentative_index);
            etaTent = etaprob(:,3);
            phiTent = phiprob(:,3);
            psNew = psprob(:,3);
            zTent = zprob(3);
            vzTent = vzprob(3);
        else
            [etaprob(:,4),phiprob(:,4),zprob(4),vzprob(4),psprob(1:(numl(tentative_index)+1),4),errortan(4,tentative_index+1)] = ...
                 getNextStep(numl(tentative_index),numl(tentative_index) + 1,dt,z(tentative_index),vz(tentative_index),etaMatPer(:, tentative_index), phiMatPer(:, tentative_index), ...
                    zs,RvTent, angleDropMP, PROBLEM_CONSTANTS, nlmaxTent); 
            [etaprob(:,2),phiprob(:,2),zprob(2),vzprob(2),psprob(1:(numl(tentative_index)-1),2),errortan(2,tentative_index+1)] = ...
                   getNextStep(numl(tentative_index),numl(tentative_index) - 1,dt,z(tentative_index),vz(tentative_index),etaMatPer(:, tentative_index), phiMatPer(:, tentative_index), ...
                        zs,RvTent, angleDropMP, PROBLEM_CONSTANTS, nlmaxTent); 
            if (abs(errortan(3, tentative_index+1)) > abs(errortan(4, tentative_index+1)) || ...
                    abs(errortan(3, tentative_index+1)) > abs(errortan(2, tentative_index+1)))
                if abs(errortan(4, tentative_index+1)) <= abs(errortan(2, tentative_index+1))
                    %Now lets check with one more point to be sure
                    [~,~,~,~,~,errortan(5,tentative_index+1)] = ...
                        getNextStep(numl(tentative_index),numl(tentative_index) + 2,dt,z(tentative_index),vz(tentative_index),etaMatPer(:, tentative_index), phiMatPer(:, tentative_index), ...
                            zs,RvTent, angleDropMP, PROBLEM_CONSTANTS, nlmaxTent); 

                    if abs(errortan(4, tentative_index+1)) < abs(errortan(5, tentative_index+1))
                        %Accept new data
                        numlTent = numl(tentative_index) + 1;
                        etaTent = etaprob(:,4);
                        phiTent = phiprob(:,4);
                        psNew = psprob(:, 4);
                        zTent = zprob(4);
                        vzTent = vzprob(4);

                    else
                        %time step too big
                        % tvec = [tvec(1:jj),tvec(jj)/2+tvec(jj+1)/2,tvec(jj+1:end)];
                        tentative_index = tentative_index-1; 
                        reduc = 1;
                    end
                else
                    %now lets check if errortan is good enough with one point
                    %less
                    [~,~,~,~,~,errortan(1,tentative_index+1)] = ...
                        getNextStep(numl(tentative_index),numl(tentative_index) - 2,dt,z(tentative_index), ...
                            vz(tentative_index),etaMatPer(:, tentative_index), phiMatPer(:, tentative_index), ...
                            zs,RvTent, angleDropMP, PROBLEM_CONSTANTS, nlmaxTent); 
                    if abs(errortan(2, tentative_index+1)) < abs(errortan(1, tentative_index+1))
                        %Accept new data
                        numlTent = numl(tentative_index) + 1;
                        etaTent = etaprob(:,2);
                        phiTent = phiprob(:,2);
                        psNew = psprob(:, 2);
                        zTent = zprob(2);
                        vzTent = vzprob(2);
                    else
                        % tvec = [tvec(1:jj),tvec(jj)/2+tvec(jj+1)/2,tvec(jj+1:end)];
                        tentative_index = tentative_index-1; 
                        reduc = 1;
                    end
                end %End of (errortan(4) < errortan(2))

            else %the same number of contact points is best    
                if errortan(3, tentative_index+1) == Inf % ALl errors are infinity
                    % tvec = [tvec(1:jj),tvec(jj)/2+tvec(jj+1)/2,tvec(jj+1:end)];
                    tentative_index = tentative_index-1; 
                    reduc = 1;
                else
                    %Accept new data
                    numlTent = numl(tentative_index);
                    etaTent = etaprob(:,3);
                    phiTent = phiprob(:,3);
                    psNew = psprob(:,3);
                    zTent = zprob(3);
                    vzTent = vzprob(3);
                end
            end %
        end

        
        if ll == 100 && reduc == 0
            % % tvec = [tvec(1:jj),tvec(jj)/2+tvec(jj+1)/2,tvec(jj+1:end)];
            tentative_index = tentative_index-1;
            reduc = 1;
        end
        
        % check whether warning was raised during loop
        [warnMsg, warnId] = lastwarn();

        % If there was a warning of dt is too small, step back one step
        if strcmp(warnId, 'MATLAB:nearlySingularMatrix') == true || dt*T_unit < 1e-10
            lastwarn('', '');

            % bc if reduc == 1, tentative index was already lowered
            if reduc == 1; tentative_index = tentative_index + 1; end
            
            t = tvec(tentative_index-1);
            %dt = dtb * 1.2;
            dt = dtb * (1 + randi([-6, 25], 1)/10); % Pseudorandom step size different from original
            time_index = 0; save_index = 0;
            tentative_index = tentative_index - 2;
            reduc = 1;

            % reinstating old variables
            ps_accepted = ps_old;
            previous_conditions{2} = previous_conditions{1};
            previous_conditions{1} = old_conditions;
            
            going_back = going_back + 1;
            disp("Warning detected. Will proceed with going back"); 
            if going_back > 100
                t = inf;
            end
        % verifying convergence of the pressure field
        elseif reduc == 0 %if there was no reduction of time step
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
            
%             err = norm([psTent;zeros(length(length(psTent)+1:nb_points),1);amplitudes_tent']-...
%                        [psNew ;zeros(length(length(psNew )+1:nb_points),1);amplitudes_new'],1);
%                               
%             if norm([psNew;amplitudes_new'],1) > 0
%                 errorP = err/norm([psNew;amplitudes_new'],1);    
%             else
%                 errorP = err;
%             end
            errorP = norm(amplitudes_tent-amplitudes_new)/norm(amplitudes_tent);
            if errorP < tolP % Finally accept solution
                
                % Advance time index (halving mechanism) and t
                time_index = time_index + 1; t = t + dt;
                if mod(time_index, 2) == 0 && resetter > 0
                    time_index = time_index / 2;
                    save_index = save_index - 1;
                    dt = 2 * dt;
                    resetter = resetter - 1;
                elseif save_index == 0 && dt ~= dtb % if mesh doesnt agree with dtb, fix it
                    if mod(t, dtb) ~= 0 
                        dt = (2-mod(t, dtb)/dtb)* dtb;
                    else
                        dt = dtb;
                    end
                end

                % Saving old variables in case we need to go back
                ps_old = ps_accepted;
                old_conditions = previous_conditions{1};
                numl(tentative_index+1) = numlTent;
                eta_accepted = etaTent;
                phi_accepted = phiTent;
                ps_accepted = psNew;
                z(tentative_index+1) = zTent;
                vz(tentative_index+1) = vzTent;
                %#---
                oscillation_amplitudes(:, tentative_index + 1) = amplitudes_new;
                etas(:, tentative_index + 1) = eta_accepted;
                Rv(tentative_index+1) = zs_from_spherical(pi, amplitudes_new);
                %amplitudes_old = amplitudes_new;
                %amplitudes_velocities_old = velocities_new;
                %B_l_ps_old = B_l_ps_new;
                
                previous_conditions{1} = previous_conditions{2};
                previous_conditions{2} = struct("deformation_amplitudes", amplitudes_new, ...
                    "deformation_velocities", velocities_new, ...
                    "dt", dt, "nb_harmonics", N, "pressure_amplitudes", B_l_ps_new, ...
                    "current_time", previous_conditions{1}.current_time + dt, ...
                    "center_of_mass", z(tentative_index+1), "center_of_mass_velocity", vz(tentative_index+1), ...
                    "nb_contact_points", numlTent);
                
                nlmax(tentative_index+1) = nlmaxTent;
                etaOri(tentative_index+1) = eta_accepted(1);
                tvec(tentative_index + 1) = t;

                % jj0 = floor(tentative_index/nsteps);
                % jj1 = round(tentative_index-jj0*nsteps);

                etaMatPer(:,tentative_index+1) = eta_accepted;
                phiMatPer(:,tentative_index+1) = phi_accepted;
                psMatPer{tentative_index+1} = ps_accepted;

                if time_index == 2^save_index
                    time_index = 0;
                    resetter = 1;
                    indexes_to_save(current_to_save) = tentative_index;
                    current_to_save = current_to_save + 1;
                end

                if zTent > 1.5 && numlTent == 0
                    tend = t;
                end
                
                %etao = eta_accepted;
                %phio = phi_accepted;
                %pso = ps_accepted;
            
                xs = dr*(0:nlmax(tentative_index+1)-1);
                zsplot = zs(1:nlmax(tentative_index+1))+RvTent+z(tentative_index+1);
                plot([-fliplr(xs(2:end)),xs],[flipud(zsplot(2:end));zsplot],'k','Linewidth',2);
                hold on
                thetaplot = linspace(0, thetaVec(end), 200);%-%-0:thetaVec(end)/400:thetaVec(end);
                
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
            dt = dt/2; time_index = 2 * time_index;
            save_index = save_index + 1;
        end %it time step was not reduced
    end % inner while   
end % Outer while

% runNumber = runNumber+1;
tstop = t;
if indexes_to_save(2) == 1; indexes_to_save = indexes_to_save(2:end); end
if t < inf
    indexes_to_save = indexes_to_save(1:(current_to_save-1));
    z = z(indexes_to_save); save('z.mat','z')
    etaOri = etaOri(indexes_to_save); save('etaOri.mat','etaOri')
    etas = etas(:, indexes_to_save); save('etas.mat', 'etas');
    etaMatPer = etaMatPer(:, indexes_to_save); save('etaMatPer.mat', etaMatPer);
    phiMatPer = phiMatPer(:, indexes_to_save); save('phiMatPer.mat','phiMatPer');
    psMatPer = psMatPer(:, indexes_to_save); save('psMatPer.mat','psMatPer');
    vz = vz(indexes_to_save); save('vz.mat','vz');
    tvec = tvec(indexes_to_save); save('tvec.mat','tvec'); 
    nlmax = nlmax(indexes_to_save); save('nlmax.mat','nlmax');
    numl = numl(indexes_to_save); save('numl.mat','numl');
    errrortan = errortan(indexes_to_save); save('errortan.mat','errortan');
    oscillation_amplitudes = oscillation_amplitudes(indexes_to_save); save('oscillation_amplitudes.mat', 'oscillation_amplitudes');
    pressure_amplitudes    = pressure_amplitudes(:, indexes_to_save);    save('pressure_amplitudes.mat', 'pressure_amplitudes');
    Rv = Rv(indexes_to_save); save('Rv.mat', 'Rv');
end