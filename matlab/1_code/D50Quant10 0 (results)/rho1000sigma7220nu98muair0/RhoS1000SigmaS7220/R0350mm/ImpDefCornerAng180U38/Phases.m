clc
% close all 
clear
load('tvec.mat')
load('Omega.mat')
load('Gamma.mat')
% load('etaOri.mat')
% load('z.mat')
% load('Ro.mat')
load('numl.mat')
flag = 0;
j = 0;%landing counter
k = 0;%take off counter
for i=1:length(numl)
    if flag == 0 && numl(i)>.5
        flag = 1;
        j = j+1;
        landing(j) = i;
    end
    if flag == 1 && numl(i)<.5
        flag = 0;
        k = k+1;
        takeOff(k) = i;
    end
end
if j>k
    landing = landing(1:end-1);
end
LandingPh = (tvec(landing)+tvec(landing-1))*80/2;
LandingPh = LandingPh - floor(LandingPh);
LandingPh = LandingPh;
figure
plot(LandingPh)
hold on
% LPh1 = mean(LandingPh(end-150:1:end))

TakeOffPh = (tvec(takeOff)+tvec(takeOff-1))*80/2;
TakeOffPh = TakeOffPh - floor(TakeOffPh);
% TakeOffPh = TakeOffPh;
TakeOffPh = TakeOffPh+.5;
TakeOffPh = TakeOffPh-floor(TakeOffPh)-.5;
plot(TakeOffPh,'k')
    
TOPh1 = mean(TakeOffPh(end-150:1:end))

% LPh2 = mean(LandingPh(end-149:2:end))
% TOPh2 = mean(TakeOffPh(end-149:2:end))

LandingPh(end-16:end)
TakeOffPh(end-16:end)

grid on

saveas(gcf,['PhasesOm0',num2str(10*Omega),'Gamma',num2str(100*Gamma),'.fig'],'fig')