clc
clear
close all

load('r.mat')
load('zs.mat','zs')
load('Int.mat')
load('tvec.mat')
load('nlmax.mat')
load('numl.mat')
load('Gamma.mat')
load('dropmass.mat')
load('g.mat')
load('dtb.mat')
load('Ro.mat','Ro')
load('vz.mat')
load('sigma.mat')
load('rho.mat')

wzero = 1;
V0 = abs(vz(1));

load('z.mat')
% load('etaOri.mat')

load('psMatPer1.mat')
psMatPer1 = psMatPer;
load('psMatPer2.mat')
psMatPer2 = psMatPer;
load('psMatPer3.mat')
psMatPer3 = psMatPer;
load('psMatPer4.mat')
psMatPer4 = psMatPer;
load('psMatPer5.mat')
psMatPer5 = psMatPer;
load('psMatPer6.mat')
psMatPer6 = psMatPer;
load('psMatPer7.mat')
psMatPer7 = psMatPer;
psMatPer = [psMatPer1,psMatPer2,psMatPer3,psMatPer4,psMatPer5,psMatPer6,psMatPer7];

numlMatPer = (psMatPer>0);
numlPer = sum(numlMatPer);
dnuml = [0,numlPer(2:end)-numlPer(1:end-1)];

load('etaMatPer1.mat')
etaMatPer1 = etaMatPer;
load('etaMatPer2.mat')
etaMatPer2 = etaMatPer;
load('etaMatPer3.mat')
etaMatPer3 = etaMatPer;
load('etaMatPer4.mat')
etaMatPer4 = etaMatPer;
load('etaMatPer5.mat')
etaMatPer5 = etaMatPer;
load('etaMatPer6.mat')
etaMatPer6 = etaMatPer;
load('etaMatPer7.mat')
etaMatPer7 = etaMatPer;
etaMatPer = [etaMatPer1,etaMatPer2,etaMatPer3,etaMatPer4,etaMatPer5,etaMatPer6,etaMatPer7];


dropWeight = dropmass*g;

ntimes = size(psMatPer,2);

tvecplot = tvec(1:size(psMatPer,2));

dtvecplot = [tvecplot(2:end)-tvecplot(1:end-1),dtb];

% zplot = z(end-ntimes:end-1);
% etaOriplot = etaOri(end-ntimes+1:end);

zmin = min(min(psMatPer));
zmax = max(max(psMatPer));

f=zeros(1,size(psMatPer,2));
fsigma = zeros(1,size(psMatPer,2));
fgrav = zeros(1,size(psMatPer,2));
for ii = 1:size(psMatPer,2)
    f(ii)=Int*psMatPer(1:nlmax,ii);
    fsigma(ii) = sum(Int(1:numlPer(ii)));
    fgrav(ii) = -rho*grav(Gamma,tvecplot(ii),g,wzero,0)*Int(1:numlPer(ii))*etaMatPer(1:numlPer(ii),ii);
end
fsigma = 2*sigma/Ro*fsigma;

fig=figure;
% G = plot(tvecplot*V0/Ro,fgrav/dropWeight,'b','LineWidth',2)
hold on
ST = plot(tvecplot*V0/Ro,fsigma/dropWeight,'--','color',[.4 .4 .4],'LineWidth',3)
% STG = plot(tvecplot*V0/Ro,(fsigma+fgrav)/dropWeight,'r','LineWidth',2)
T = plot(tvecplot*V0/Ro,f/dropWeight,'k','LineWidth',4)
% Leg = legend([G STG T],'Grav.','Surf. Tension+Grav.','Total')
% set(Leg,'FontSize',16,'Location','NorthEast')


% plot(tvecplot*V0/Ro,numl(1:length(tvecplot)),'color',[.5 .5 .5],'LineWidth',2)

Integral = f*dtvecplot'
% ref = dropmass*g*1/40

hold on
grid on



% for ii=1:length(dnuml)
%     if dnuml(ii) > .5
%         plot(80*[tvecplot(ii) tvecplot(ii)],[-1 10],'color',[.5 .5 .5])
%     elseif dnuml(ii) < -.5
%         plot(80*[tvecplot(ii) tvecplot(ii)],[-1 10],'--','color',[.5 .5 .5])
%     end
% end

xlabel('   $tV_0/R_o$   ','interpreter','LaTeX','FontSize',24)
ylabel('   $\frac{F_{\uparrow}}{mg}\ \ \ $    ','interpreter','LaTeX','FontSize',24,'Rotation',0)
% title(['a) $\Omega$ = ',num2str(Omega),', $\Gamma$ = ',num2str(Gamma)],'interpreter','LaTeX','FontSize',24)
set(gca,'xlim',[0 22],'ylim',[0 70],'FontName','Times','FontSize',24)%,'Xtick',[-.5:.25:.5],'Ytick',[]);
title('   b   ','FontName','times','FontSize',24)
box on
saveas(gcf,'forces21ExpLee.fig','fig')
print(fig,'-depsc','-r300',['Forces21ExpLee.eps'])



% figure
% for ii = size(etaMatPer,2)/2:15:size(etaMatPer,2)
%     plot([fliplr(-1*r),r(2:end)],(ii-1)*delta+[flipud(etaMatPer(:,ii));etaMatPer(2:end,ii)],'k','LineWidth',2)
%     grid on
%     hold on
% %     plot(r,zplot(ii)+zs)
% %     axis equal
%     set(gca,'xlim',[-.5 .5],'ylim',[.08 .19]);
% end






