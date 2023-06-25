    U0 = 35; save('U0.mat','U0')%impact velocity in cm/s
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
    load('angleDrop.mat','angleDrop')
    load('Int.mat','Int')
    load('DTNnew345nr2500D100refp10.mat','DTNnew345')
    DTN = DTNnew345;
    clear DTNnew345
    load('zs.mat','zs')
    cd ..
    load('nu.mat','nu')
    load('g.mat','g')
    load('sigma.mat','sigma')
    load('rho.mat','rho')
    
    cd(['NewR0',num2str(Ro*10000),'mm'])
    cd(['Rho',num2str(rhoS*1000)])
    cd(['U',num2str(U0)])

    tiempoComp = zeros(1,10); %just to check how long it takes to solve the first ten periods
    
    %Unit of time
    T = Ro/U0; save('T.mat','T')%base time is seconds

    %Dimensionless numbers that depend on U0
    Re = Ro*U0/nu; save('Re.mat','Re')
    Fr = U0^2/(g*Ro); save('Fr.mat','Fr')
    We = rho*U0^2*Ro/sigma; save('We.mat','We')
    WeS = rhoS*Ro*U0^2/sigma;save('WeS.mat','WeS')
    Cang = 155*pi/180; save('Cang.mat','Cang')
    
    %Physical parameters
    tend = 15; save('tend.mat','tend')%Earliest possible end of simulation in characteristic units
    
    %Numerical Simulation parameters
    nsteps = 40; save('nsteps.mat','nsteps')%number of timesteps in one unit of time
    dtb = 1/nsteps; save('dtb.mat','dtb')%basic timestep
    
    %Inintial conditions for the fluid
    t = 0;
    steps = ceil((tend-t)/dtb); %estimated total number of timesteps
    etao = zeros(nr,1); %initial surface elevation
    phio = zeros(nr,1); %initial surface potential
    pso = zeros(nlmax,1); %initial pressure distribution on surface
    
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
    
    %Initial conditions for the drop
    z(1) = 1/20; %in dimensionless units
    DelV = (1-sqrt(1-6*z(1)/Fr))/3;
    vz(1) = -1+DelV; %in dimesionless units
    
    j = 0;%iteration counter

    errortan = zeros(5,steps+1); save('errortan.mat','errortan')

    eta = etao;
    phi = phio;
    ps = pso;
    
    eta1 = eta;
    phi1 = phi;
    ps1 = ps;
    
    etaMatPer = zeros(length(eta),nsteps);
    phiMatPer = zeros(length(phi),nsteps);
    psMatPer = zeros(nlmax,nsteps);
    
    etaMatPer(:,1) = eta;
    phiMatPer(:,1) = phi;
    psMatPer(:,1) = ps;
    
    j1 = 1; %partial results savings  counter