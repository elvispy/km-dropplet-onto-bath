clear
close all
clc

tstart = tic;
%data in cgs
tmax = 15000;

load('runNumber.mat','runNumber')

% if runNumber == 0
    U0 = 40; save('U0.mat','U0')%impact velocity in cm/s
    cd ..
    load('rhoS.mat','rhoS')
    load('dropmass.mat','dropmass')
    load('Ma.mat','Ma')
    cd ..
    load('Ro.mat','Ro')
    load('alpha.mat','alpha')
    load('nr.mat','nr')
    load('nlmax.mat','nlmax')
    load('dr.mat','dr')
    load('Delta.mat','Delta')
    load('angleDropMP.mat','angleDropMP')
    load('IntMat.mat','IntMat')
    load('DTNnew345nr2500D100refp10.mat','DTNnew345')
    DTN = DTNnew345;
    clear DTNnew345
    load('zs.mat','zs')
    load('xdrop.mat')
    load('xplot.mat')
    load('zdrop.mat')
    cd ..
    load('nu.mat','nu')
    load('g.mat','g')
    load('sigma.mat','sigma')
    load('rho.mat','rho')
    
    cd(['CuspR0',num2str(Ro*10000),'mm'])
    cd(['Rho',num2str(rhoS*1000)])
    cd(['ImpDefAng0U',num2str(U0)])

    tiempoComp = zeros(1,10); %jjust to check how long it takes to solve the first ten periods
    
    %Unit of time
    T = Ro/U0; save('T.mat','T')%base time is seconds

    %Dimensionless numbers that depend on U0
    Re = Ro*U0/nu; save('Re.mat','Re')
    Fr = U0^2/(g*Ro); save('Fr.mat','Fr')
    We = rho*U0^2*Ro/sigma; save('We.mat','We')
    WeS = rhoS*Ro*U0^2/sigma;save('WeS.mat','WeS')
    Cang = pi; save('Cang.mat','Cang')
    
    %Physical parameters
    tend = 40; save('tend.mat','tend')%Earliest possible end of simulation in characteristic units
    
    %Numerical Simulation parameters
    nsteps = 40; save('nsteps.mat','nsteps')%number of timesteps in one unit of time
    dtb = 1/nsteps; save('dtb.mat','dtb')%basic timestep
%     dt = dtb; %seting the current timestep to the basic value
    
    %Inintial conditions for the fluid
    t = 0;
    steps = ceil((tend-t)/dtb); %estimated total number of timesteps
    etao = zeros(nr,1); %initial surface elevation
    phio = zeros(nr,1); %initial surface potential
%     pso = zeros(nlmax,1); %initial pressure distribution on surface
    
    %Zeroing result storing variables
    etaOri = zeros(1,steps+1);
    save('etaOri.mat','etaOri')
    z = zeros(1,steps+1);
    save('z.mat','z')
    vz = zeros(1,steps+1);
    save('vz.mat','vz')
    numl = zeros(1,steps+1);
    save('numl.mat','numl')
    tvec = t:dtb:tend+1; %giving extra time just in case the simulation needs to run longer
    save('tvec.mat','tvec')
    Rv = zeros(1,steps+1);
    V = zeros(1,steps+1);
    nlmax = zeros(1,steps+1);
    
    %Drop deformation parameters
    alpha = .01/Re;%3.8/Re;
    beta = .01/We;%5.84/We;
    tolP = .01; %error tolerance for the pressure field
    
    %Initial conditions for the drop
    Rv(1) = 1;
    V(1) = 0;
    
    %Consecuences of the initial conditions (assuming ellipse)
    Rh(1) = sqrt(1/Rv(1));
    nlmax(1) = floor(Rh(1)/dr)+1;
    zs(1:nlmax(1)) = Rv(1)-Rv(1)*sqrt(1-Rv(1)*(0:dr:(nlmax(1)-1)*dr).^2);
    zs(nlmax(1)+1:end) = 1000;
    if floor(Rh/dr) == Rh/dr
        tanDrop(1:nlmax(1)-1) = Rv(1)^2*(.5*dr:dr:(nlmax(1)-1)*dr)./...
            sqrt(1-Rv(1)*(.5*dr:dr:(nlmax(1)-1)*dr).^2);
        angleDropMP(1:nlmax(1)-1) = atan(tanDrop(1:nlmax(1)-1));
        angleDropMP(nlmax(1)) = pi/2;%#ok
    else
        tanDrop(1:nlmax(1)) = Rv(1)^2*(.5*dr:dr:(nlmax(1))*dr)./sqrt(1-Rv(1)*(.5*dr:dr:(nlmax(1))*dr).^2);
        angleDropMP(1:nlmax(1)) = atan(tanDrop(1:nlmax(1)));
    end
    z(1) = Rv(1); %in dimensionless units
    vz(1) = -1; %in dimesionless units
    
    jj = 0;%iteration counter

    errortan = zeros(5,steps+1); save('errortan.mat','errortan')

    eta = etao;
    phi = phio;
    
    eta1 = eta;
    phi1 = phi;
    ps1 = [];
    
    etaMatPer = zeros(length(eta),nsteps);
    phiMatPer = zeros(length(phi),nsteps);
    psMatPer = cell(1,nsteps);
    
    etaMatPer(:,1) = eta;
    phiMatPer(:,1) = phi;
    psMatPer{1} = [];
    
    jj1 = 1; %partial results savings  counter

% elseif runNumber>0
%     load('U0.mat','U0')%impact velocity in cm/s
%     cd ..
%     load('rhoS.mat','rhoS')
%     load('dropmass.mat','dropmass')
%     load('Ma.mat','Ma')
%     cd ..
%     load('Ro.mat','Ro')
%     load('alpha.mat','alpha')
%     load('nr.mat','nr')
%     load('nlmax.mat','nlmax')
%     load('dr.mat','dr')
%     load('Delta.mat','Delta')
%     load('angleDropMP.mat','angleDropMP')
%     load('IntMat.mat','IntMat')
%     load('DTNnew345nr2500D100refp10.mat','DTNnew345')
%     DTN = DTNnew345;
%     clear DTNnew345
%     load('zs.mat','zs')
%     cd ..
%     load('nu.mat','nu')
%     load('g.mat','g')
%     load('sigma.mat','sigma')
%     load('rho.mat','rho')
%     cd(['NewR0',num2str(Ro*10000),'mm'])
%     cd(['Rho',num2str(rhoS*1000)])
%     cd(['U',num2str(U0)])
%     
%     tiempoComp = zeros(1,10); %just to check how long it takes to solve the first ten periods
%     
%     %Unit of time
%     load('T.mat','T')%base time is seconds
%     %Dimensionless numbers that depend on U0
%     load('Re.mat','Re')
%     load('Fr.mat','Fr')
%     load('We.mat','We')
%     load('WeS.mat','WeS')
%     load('Cang.mat','Cang')
%     
%     %Physical parameters
%     load('tend.mat','tend')%Earliest possible end of simulation in characteristic units
%     
%     %Numerical Simulation parameters
%     load('nsteps.mat','nsteps')%number of timesteps in one unit of time
%     load('dtb.mat','dtb')%basic timestep
%     
%     %Zeroing result storing variables
%     load(['etao',num2str(runNumber),'.mat'],'etao')
%     load(['phio',num2str(runNumber),'.mat'],'phio')
%     load(['pso',num2str(runNumber),'.mat'],'pso')
%     load('etaOri.mat','etaOri')
%     load('z.mat','z')
%     load('vz.mat')
%     load('numl.mat')
%     load('tvec.mat')
% 
%     %Inintial conditions
%     load(['tstop',num2str(runNumber),'.mat'],'tstop') 
%     t = tstop;
%     load(['jjstop',num2str(runNumber),'.mat'],'jjstop')
%     jj = jjstop;
%     
%     eta = etao;
%     phi = phio;
%     ps = pso;
%     
%     etaMatPer = zeros(length(eta),nsteps);
%     phiMatPer = zeros(length(phi),nsteps);
%     psMatPer = cell(1,nsteps);
%     
%     etaMatPer(:,1) = eta;
%     phiMatPer(:,1) = phi;
%     psMatPer{1} = ps;
% 
%     %Drop geometry
%     load('xdrop.mat','xdrop')
%     load('zdrop.mat','zdrop')
%     load('angleDropMP.mat','angleDropMP')
%     load('errortan.mat','errortan')
%     load('IntMat.mat','IntMat')
%     
%     jj1 = 1;
% 
% elseif runNumber < 0
%     exit
% end

%Main Loop
while t<tend || jj1>.5
    jj = jj+1;
    t = tvec(jj+1);
    dt = t - tvec(jj);

    etaprob = zeros(nr,5);
    phiprob = zeros(nr,5);
    vzprob = zeros(1,5);
    zprob = zeros(1,5);
    errortan(:,jj+1) = 4*ones(5,1);
    ps0 = ps1;
    psTent = ps0; %Tentative pressure distribution
    A = [1+alpha*dt/2 dt/2; -beta*dt/2 1]/(beta*dt^2/4+(1+alpha*dt/2));%lhs ode matrix inverted
    %rhs ode matrix
    B = [1+alpha*dt/2-beta*dt^2/4 dt; -beta*dt beta*dt^2/4+(1-alpha*dt/2)]/(beta*dt^2/4+(1+alpha*dt/2));
    if numl(jj) == 0%integrating pressure
        Force = 0;
    else
        Force = IntMat(numl(jj),1:length(ps0))*ps0;
    end
    TentVec = B*[Rv(jj); V(jj)]+A*[0; beta*dt-3*dt*Force/(4*pi)];%tentative Rv and dt Rv
    RvTent = TentVec(1);
    VTent = TentVec(2);%tentative velocity of radius deformation
    RhTent = sqrt(1/RvTent);%tentative radius deformation
    nlmaxTent = floor(RhTent/dr)+1;%Tentative max number of mesh points covered by the deformed droplet
    %equation for the bottom of the tentative droplet shape
    zs(1:nlmaxTent) = RvTent-RvTent*sqrt(1-RvTent*(0:dr:(nlmaxTent-1)*dr).^2);
    zs(nlmaxTent+1:end) = 1000;
    %tangent direction to the tentative shape
    if floor(RhTent/dr) == RhTent/dr
        tanDrop(1:nlmaxTent-1) = RvTent^2*(.5*dr:dr:(nlmaxTent-1)*dr)./...
                                sqrt(1-RvTent*(.5*dr:dr:(nlmaxTent-1)*dr).^2);
        angleDropMP(1:nlmaxTent-1) = atan(tanDrop(1:nlmaxTent-1));
        angleDropMP(nlmaxTent) = pi/2;%#ok
    else
        tanDrop(1:nlmaxTent) = RvTent^2*(.5*dr:dr:(nlmaxTent)*dr)./...
                                sqrt(1-RvTent*(.5*dr:dr:(nlmaxTent)*dr).^2);
        angleDropMP(1:nlmaxTent) = atan(tanDrop(1:nlmaxTent));
    end
    
    psprob = zeros(nlmaxTent,5);%zeroing the vector of potential pressures
    
    errorP = 1; %tolerance for force discrepancy
    reduc = 0; %indicator of whether there was a reduction in the time-step size or not
    
    while abs(errorP)>=tolP && reduc == 0
        if numl(jj) < .5 %i.e. if previously in flight
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
                    We,Ma,zs,IntMat(1,:),angleDropMP,Cang,WeS,RvTent);
                co = find(numl(jj:-1:1)~=2,1);
                [~,~,~,~,~,errortan(5,jj+1)] = ...    
                    solvenDDCusp(numl(jj-co+1),2,dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(2,:),angleDropMP,Cang,WeS,RvTent);
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
                    We,Ma,zs,IntMat(1,:),angleDropMP,Cang,WeS,RvTent);
                co = find(numl(jj:-1:1)~=2,1);
                [etaprob(:,4),phiprob(:,4),zprob(4),vzprob(4),psprob(1:2,4),errortan(4,jj+1)] = ...    
                    solvenDDCusp(numl(jj-co+1),2,dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(2,:),angleDropMP,Cang,WeS,RvTent);
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
                        We,Ma,zs,IntMat(3,:),angleDropMP,Cang,WeS,RvTent);
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
                    We,Ma,zs,IntMat(2,:),angleDropMP,Cang,WeS,RvTent);
                co = find(numl(jj:-1:1)~=1,1);
                [etaprob(:,2),phiprob(:,2),zprob(2),vzprob(2),psprob(1,2),errortan(2,jj+1)] = ...    
                    solvenDDCusp(numl(jj-co+1),1,dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(1,:),angleDropMP,Cang,WeS,RvTent);
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
                        We,Ma,zs,IntMat(3,:),angleDropMP,Cang,WeS,RvTent);
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
                            We,Ma,zs,IntMat(4,:),angleDropMP,Cang,WeS,RvTent);
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
        elseif numl(jj)>2.5 && numl(jj)<nlmax(jj)-1.5 %if the last number of contact points was far from the boundaries
            co = find(numl(jj:-1:1)~=numl(jj),1);
            [etaprob(:,3),phiprob(:,3),zprob(3),vzprob(3),psprob(1:numl(jj),3),errortan(3,jj+1)] = ...    
                solvenDDCusp(numl(jj-co+1),numl(jj),dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                We,Ma,zs,IntMat(numl(jj),:),angleDropMP,Cang,WeS,RvTent);
            co = find(numl(jj:-1:1)~=numl(jj)-1,1);
            [etaprob(:,2),phiprob(:,2),zprob(2),vzprob(2),psprob(1:numl(jj)-1,2),errortan(2,jj+1)] = ...    
                solvenDDCusp(numl(jj-co+1),numl(jj)-1,dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                We,Ma,zs,IntMat(numl(jj)-1,:),angleDropMP,Cang,WeS,RvTent);
            if abs(errortan(2,jj+1)) < abs(errortan(3,jj+1))
                co = find(numl(jj:-1:1)~=numl(jj)-2,1);
                [~,~,~,~,~,errortan(1,jj+1)] = ...
                    solvenDDCusp(numl(jj-co+1),numl(jj)-2,dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(numl(jj)-2,:),angleDropMP,Cang,WeS,RvTent);
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
                    We,Ma,zs,IntMat(numl(jj)+1,:),angleDropMP,Cang,WeS,RvTent);
                if abs(errortan(3,jj+1))<abs(errortan(4,jj+1))
                    numlTent = numl(jj);
                    etaTent = etaprob(:,3);
                    phiTent = phiprob(:,3);
                    psNew = psprob(:,3);
                    zTent = zprob(3);
                    vzTent = vzprob(3);
                else
                    co = find(numl(jj:-1:1)~=numl(jj)+2,1);
                    [~,~,~,~,~,errortan(5,jj+1)] = ...
                        solvenDDCusp(numl(jj-co+1),numl(jj)+2,dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                        We,Ma,zs,IntMat(numl(jj)+2,:),angleDropMP,Cang,WeS,RvTent);
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
        elseif numl(jj) > nlmax-1.5 && numl(jj) < nlmax(jj)-.5 %i.e. if last number of contacted points is nlmax-1
            co = find(numl(jj:-1:1)~=numl(jj),1);
            [etaprob(:,3),phiprob(:,3),zprob(3),vzprob(3),psprob(1:numl(jj),3),errortan(3,jj+1)] = ...    
                solvenDDCusp(numl(jj-co+1),numl(jj),dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                We,Ma,zs,IntMat(numl(jj),:),angleDropMP,Cang,WeS,RvTent);
            co = find(numl(jj:-1:1)~=numl(jj)-1,1);
            [etaprob(:,2),phiprob(:,2),zprob(2),vzprob(2),psprob(1:numl(jj)-1,2),errortan(2,jj+1)] = ...    
                solvenDDCusp(numl(jj-co+1),numl(jj)-1,dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                We,Ma,zs,IntMat(numl(jj)-1,:),angleDropMP,Cang,WeS,RvTent);
            if abs(errortan(2,jj+1))<abs(errortan(3,jj+1))
                co = find(numl(jj:-1:1)~=numl(jj)-2,1);
                [~,~,~,~,~,errortan(1,jj+1)] = ...    
                    solvenDDCusp(numl(jj-co+1),numl(jj)-2,dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(numl(jj)-2,:),angleDropMP,Cang,WeS,RvTent);
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
                    We,Ma,zs,IntMat(numl(jj)+1,:),angleDropMP,Cang,WeS,RvTent);
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
        else %i.e. if last number of contact points was nlmax
            co = find(numl(jj:-1:1)~=numl(jj),1);
            [etaprob(:,3),phiprob(:,3),zprob(3),vzprob(3),psprob(1:numl(jj),3),errortan(3,jj+1)] = ...    
                solvenDDCusp(numl(jj-co+1),numl(jj),dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                We,Ma,zs,IntMat(numl(jj),:),angleDropMP,Cang,WeS,RvTent);
            co = find(numl(jj:-1:1)~=numl(jj)-1,1);
            [etaprob(:,2),phiprob(:,2),zprob(2),vzprob(2),psprob(1:numl(jj)-1,2),errortan(2,jj+1)] = ...    
                solvenDDCusp(numl(jj-co+1),numl(jj)-1,dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                We,Ma,zs,IntMat(numl(jj)-1,:),angleDropMP,Cang,WeS,RvTent);
            if abs(errortan(2,jj+1)) < abs(errortan(3,jj+1))
                co = find(numl(jj:-1:1)~=numl(jj)-2,1);
                [~,~,~,~,~,errortan(1,jj+1)] = ...
                    solvenDDCusp(numl(jj-co+1),numl(jj)-2,dt,z(jj),vz(jj),etao,phio,nr,dr,Re,Delta,DTN,Fr,...
                    We,Ma,zs,IntMat(numl(jj)-2,:),angleDropMP,Cang,WeS,RvTent);
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
        end
        

        % verifying convergence of the pressure field
        if reduc == 0
            nume = max(length(psTent),length(psNew));
            errP = norm([psTent;zeros(length(length(psTent)+1:nume),1)]-...
                          [psNew ;zeros(length(length(psNew )+1:nume),1)],1);
            if errP > 0 %I suspect that the argument of this if should be the quantity that errP is divided by
                errorP = errP/norm([psNew;zeros(length(length(psNew)+1:nume),1)],1);
            else
                errorP = errP;
            end
            if errorP < tolP
                numl(jj+1) = numlTent;
                eta1 = etaTent;
                phi1 = phiTent;
                ps1 = psNew;
                z(jj+1) = zTent;
                vz(jj+1) = vzTent;
                Rv(jj+1) = RvTent;
                V(jj+1) = VTent;
                nlmax(jj+1) = nlmaxTent;
                etaOri(jj+1) = eta1(1);


                eta = eta1;
                phi = phi1;

                jj0 = floor(jj/nsteps);
                jj1 = round(jj-jj0*nsteps);
            %     if jj<10*nsteps || jj>=222*nsteps
                etaMatPer(:,jj1+1) = eta;
                phiMatPer(:,jj1+1) = phi;
                psMatPer{jj1+1} = ps1;
                if jj1 == nsteps-1
                    if runNumber == 0
                        tiempoComp(jj0+1)=toc(tstart);
                    end
                    save(['etaMatPer',num2str(jj0+1),'.mat'],'etaMatPer')
                    save(['phiMatPer',num2str(jj0+1),'.mat'],'phiMatPer')
                    save(['psMatPer',num2str(jj0+1),'.mat'],'psMatPer')
            %         save('qq.mat','qq')
                    save('etaOri.mat','etaOri')
                    save('Rv.mat','Rv')
                    save('V.mat')
                    save('z.mat','z')
                    save('vz.mat','vz')
                    save('tvec.mat','tvec')
                    save('numl.mat','numl')
                    save('errortan.mat','errortan')
                end
            %     end
                etao = eta1;
                phio = phi1;
            %     pso = ps1;

                etaplot=[flipud(eta(2:nr));eta];
                xs = [xplot(nr-nlmax(jj+1)+2:nr+nlmax(jj+1)),...
                    fliplr(xplot(nr-nlmax(jj+1)+2:nr+nlmax(jj+1))),...
                    xplot(nr-nlmax(jj+1)+2)];
                zsplot = [(z(jj+1)-Rv(jj+1))+[flipud(zs(2:nlmax(jj+1)));zs(1:nlmax(jj+1))];...
                            flipud((z(jj+1)+Rv(jj+1))-[flipud(zs(2:nlmax(jj+1)));...
                            zs(1:nlmax(jj+1))]);(z(jj+1)-Rv(jj+1))+zs(nlmax(jj+1))];
                plot(xs,zsplot,'k','Linewidth',4);
                hold on
                width = 200;
                plot(xplot(nr+1-width:nr+1+width),etaplot(nr-width:nr+width),'LineWidth',4);
                hold off
                axis equal
                title(['   t = ',num2str(t),'   ','nl = ',num2str(numl(jj+1))],'FontSize',16);
                grid on
                drawnow;
            else
                psTent = psNew;
                if numlTent > 0
                    TentVec = B*[Rv(jj); V(jj)]+...
                                A*[0; beta*dt-3*dt*(Force+IntMat(numlTent,1:length(psTent))*psTent)/(8*pi)];
                else
                    TentVec = B*[Rv(jj); V(jj)]+...
                                A*[0; beta*dt-3*dt*(Force)/(8*pi)];
                end
                RvTent = TentVec(1);
                VTent = TentVec(2);
                RhTent = sqrt(1/RvTent);
                nlmaxTent = floor(RhTent/dr)+1;
                zs(1:nlmaxTent) = RvTent-RvTent*sqrt(1-RvTent*(0:dr:(nlmaxTent-1)*dr).^2);
                zs(nlmaxTent+1:end) = 1000;
                if floor(RhTent/dr) == RhTent/dr
                    tanDrop(1:nlmaxTent-1) = RvTent^2*(.5*dr:dr:(nlmaxTent-1)*dr)./...
                                            sqrt(1-RvTent*(.5*dr:dr:(nlmaxTent-1)*dr).^2);
                    angleDropMP(1:nlmaxTent-1) = atan(tanDrop(1:nlmaxTent-1));
                    angleDropMP(nlmaxTent) = pi/2;%#ok
                else
                    tanDrop(1:nlmaxTent) = RvTent^2*(.5*dr:dr:(nlmaxTent)*dr)./...
                                            sqrt(1-RvTent*(.5*dr:dr:(nlmaxTent)*dr).^2);
                    angleDropMP(1:nlmaxTent) = atan(tanDrop(1:nlmaxTent));
                end
            end
        end
    end
    

    tComp = toc(tstart);
    if jj1==0 && tComp > tmax
        runNumber = runNumber+1;
        save('runNumber.mat','runNumber')
        tstop = t;
        save(['tstop',num2str(runNumber),'.mat'],'tstop')
        jjstop = jj;
        save(['jjstop',num2str(runNumber),'.mat'],'jjstop')
        save(['etao',num2str(runNumber),'.mat'],'etao')
        save(['phio',num2str(runNumber),'.mat'],'phio')
%         save(['pso',num2str(runNumber),'.mat'],'pso')

        zrestart = z(jj+1);
        vzrestart = vz(jj+1);
        trestart = tvec(jj+1);
        numlrestart = numl(jj+1);
        save(['zrestart',num2str(runNumber),'.mat'],'zrestart')
        save(['vzrestart',num2str(runNumber),'.mat'],'vzrestart')
        save(['trestart',num2str(runNumber),'.mat'],'trestart')
        save(['numlrestart',num2str(runNumber),'.mat'],'numlrestart')

        save('etaOri.mat','etaOri')
        save('Rv.mat','Rv')
        save('V.mat')
        save('z.mat','z')
        save('vz.mat','vz')
        save('tvec.mat','tvec')
        save('numl.mat','numl')
        save('errortan.mat','errortan')
        if runNumber == 1
            save('tiempoComp.mat','tiempoComp')
        end
        exit
    end
end

runNumber = runNumber+1;
tstop = t;
save(['tstop',num2str(runNumber),'.mat'],'tstop')
jjstop = jj;
save(['jjstop',num2str(runNumber),'.mat'],'jjstop')
save(['etao',num2str(runNumber),'.mat'],'etao')
save(['phio',num2str(runNumber),'.mat'],'phio')
% save(['pso',num2str(runNumber),'.mat'],'pso')

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
save('Rv.mat')
save('V.mat')
save('vz.mat','vz')
save('tvec.mat','tvec')
save('numl.mat','numl')
save('errortan.mat','errortan')

% runNumber = -10;
% save('runNumber.mat','runNumber')
% exit






