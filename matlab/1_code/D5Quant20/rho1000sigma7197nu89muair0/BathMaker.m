clear
close all
clc

%Physical properties of the bath's fluid
rho = 1; save('rho.mat','rho')%Density of fluid bath in cgs
g = 980; save('g.mat','g')%value of gravity in cgs
sigma = 71.97; save('sigma.mat','sigma')%Surface tension in cgs
nu = 8.94E-3; save('nu.mat','nu')%kinematic viscosity of fluid in cgs (8.94E-3)

%Physical properties of the air
muair = 0; save('muair.mat','muair')%Dynamic viscosity of air in cgs (typically 1.84*10^(-4), here neglected)

