clear
close all
clc

rhoS = 1.0; save('rhoS.mat','rhoS')%ball density in gr/cm^3
sigmaS = 71.97; save('sigmaS.mat','sigmaS')

cd ..
load('rho.mat','rho')

cd(['RhoS',num2str(1000*rhoS),'SigmaS',num2str(sigmaS)])
Ma = 4*pi*rhoS/(3*rho); save('Ma.mat','Ma') %Dimensionless mass of the droplet
Ra = rhoS/rho; save('Ra.mat','Ra') %Density ratio of the two fluids
