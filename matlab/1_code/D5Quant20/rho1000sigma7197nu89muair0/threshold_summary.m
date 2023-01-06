close all
clear
clc

rho = 1;
sigma = 71.97;
g = 980;

%Ro = 0.083 cm
Bo_Ro_083 = rho*g*.083^2/sigma;

densidades_Ro_083 = [.25 .5 .8 1 1.2];
D_Ro_083 = densidades_Ro_083;

Velocities_Ro_083 = [7.5 12.5 16.5 18.5 19.5];
We_Ro_083 =rho*.083*Velocities_Ro_083.^2/sigma;

figure
plot(log(D_Ro_083),log(We_Ro_083),'-x','LineWidth',2,'MarkerSize',10)
grid on
hold on
coeffs = polyfit(log(D_Ro_083),log(We_Ro_083),1)
plot(log(D_Ro_083),coeffs(1)*log(D_Ro_083)+coeffs(2),'--+')

A1 = exp(coeffs(2))/(Bo_Ro_083^1.79)


%rhoS = 1.2

