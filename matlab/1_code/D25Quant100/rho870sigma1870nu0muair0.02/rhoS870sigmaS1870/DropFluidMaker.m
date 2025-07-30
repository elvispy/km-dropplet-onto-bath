
rhoS = 0.87; save('rhoS.mat','rhoS')%ball density in gr/cm^3
sigmaS = 18.7; save('sigmaS.mat','sigmaS');

cd ..
load('rho.mat','rho');

cd(sprintf("Rhos%gSigmaS%g", 1000*rhoS, 100*sigmaS));
Ma = 4*pi*rhoS/(3*rho); save('Ma.mat','Ma'); %Dimensionless mass of the droplet
Ra = rhoS/rho; save('Ra.mat','Ra'); %Density ratio of the two fluids
