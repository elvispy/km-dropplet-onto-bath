
%Physical properties of the bath's fluid
rho = 1; save('rho.mat','rho')%Density of fluid bath in cgs
g = 981; save('g.mat','g')%value of gravity in cgs
sigma = 72.2; save('sigma.mat','sigma')%Surface tension in cgs. Updated to 72.22 dyne/cm to match Alventosa et al (2022)
nu = 9.78E-3; save('nu.mat','nu')%Updated to 9.78E-3 cm2/s to match Alventosa et al (2022)

%Physical properties of the air
muair = 0; save('muair.mat','muair')%Dynamic viscosity of air in cgs (typically 1.84*10^(-4), here neglected)

