clear
clc
close all

load('tvec.mat','tvec')
load('U0.mat','U0')
load('numl.mat','numl')
cd ..
load('rhoS.mat')
cd ..
load('Ro.mat','Ro')
load('dr.mat')
load('IntMat.mat')
cd(['Rho',num2str(rhoS*1000)])
load('Ma.mat')
cd(['ImpDefAng0U',num2str(U0)])
load('Fr.mat')
load('z.mat','z')
load('vz.mat')
load('etaOri.mat','etaOri')
load('Rv.mat')
load('dtb.mat')

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
load('psMatPer8.mat')
psMatPer8 = psMatPer;
load('psMatPer9.mat')
psMatPer9 = psMatPer;
load('psMatPer10.mat')
psMatPer10 = psMatPer;
load('psMatPer11.mat')
psMatPer11 = psMatPer;
load('psMatPer12.mat')
psMatPer12 = psMatPer;
load('psMatPer13.mat')
psMatPer13 = psMatPer;
load('psMatPer14.mat')
psMatPer14 = psMatPer;
load('psMatPer15.mat')
psMatPer15 = psMatPer;
load('psMatPer16.mat')
psMatPer16 = psMatPer;
psMatPer = [psMatPer1,psMatPer2,psMatPer3,psMatPer4,psMatPer5,psMatPer6,...
            psMatPer7,psMatPer8,psMatPer9,psMatPer10,psMatPer11,...
            psMatPer12,psMatPer13,psMatPer14,psMatPer15,psMatPer16];

fi = figure;

% tExpDan = (data(:,1)-.007475)*Vo/Ro;
% zExpDan = (data(:,2)+0.6403)/(10*Ro);
% devExpDan = (data(:,3))/(10*Ro);
% Exp = fill([tExpDan;flipud(tExpDan)],[zExpDan+devExpDan;flipud(zExpDan-devExpDan)],[.5 .5 .5],'EdgeColor',[.5 .5 .5])
hold on
% ExpLine = plot(tExpDan,zExpDan,'k')

FreeSurf = plot(tvec(1:length(etaOri)),etaOri,'color',[0 0 1],'LineWidth',4)

South = plot(tvec(1:length(z)-1),z(1:end-1)-Rv,'color',[ 0.4660    0.6740    0.1880],'LineWidth',4)

Center = plot(tvec(1:length(z)),z,'k','LineWidth',4)
set(gca,'FontSize',16,'xlim',[0 16],'ylim',[-2 8])
xlabel('   $tV_0/R_o $   ','interpreter','LaTeX','FontSize',24)
ylabel('   $\frac{z}{R_o}\ \ \ $    ','interpreter','LaTeX','FontSize',32,'Rotation',0)
PressedRad = plot(tvec(1:length(z)),dr*numl)
% Force = 
title(['$R_0\, =\ $',num2str(round(10000*Ro)/1000),'$mm$, $\rho_s\, =\ $',num2str(rhoS),'$gr/cm^3$, $V_0\, =\ $',num2str(U0),'$cm/s$'],'interpreter','LaTeX','FontSize',20)
dropWeight = Ma/Fr;

ntimes = size(psMatPer,2);

tvecplot = tvec(1:size(psMatPer,2));

dtvecplot = [tvecplot(2:end)-tvecplot(1:end-1),dtb];


f=zeros(1,size(psMatPer,2));

% fig=figure;
for ii = 1:size(psMatPer,2)
    if numl(ii) == 0
        f(ii) = 0;
    else
        f(ii) = IntMat(numl(ii),1:numl(ii))*psMatPer{ii}(1:numl(ii));
    end
end
plot(tvecplot,f/dropWeight,'LineWidth',2)
hold on
plot(tvecplot,dr*numl(1:length(tvecplot)),'LineWidth',2)

Integral = f*dtvecplot'
% daspect([.83 .73/5 1])
grid on
% %
saveas(gcf,['CenterLineRadius',num2str(10*Ro),'mm.fig'],'fig')
print(fi,'-depsc','-r300',['CenterLineRadius',num2str(10*Ro),'mm.eps'])

index1 = find(z<0,1);
index2 = find(z(index1:end)>0,1);

tImpact = (tvec(index1)+tvec(index1-1))/2; save('tImpact.mat','tImpact')
Uo = (vz(index1)+vz(index1-1))/2;save('Uo.mat','Uo')
tend = (tvec(index1+index2-1)+tvec(index1+index2-2))/2;save('tend.mat','tend')
Uend = (vz(index1+index2-1)+vz(index1+index2-2))/2;save('Uend.mat','Uend')
tcont = tend-tImpact;save('tcont.mat','tcont')
CRref = -Uend/Uo; save('CRref.mat','CRref')