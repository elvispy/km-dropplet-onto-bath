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

ntimes = size(psMatPer,2);

tvecplot = tvec(1:size(psMatPer,2));

dtvecplot = [tvecplot(2:end)-tvecplot(1:end-1),dtb];

zmin = -2;
zmax = 100;

% f=zeros(1,size(psMatPer,2));

fig=figure;
for ii = 1:size(psMatPer,2)
    nl = numl(ii);
    nlmax = 2*size(IntMat,1);
    plot(dr*[0:nl-1,nl-.5,nl:nlmax-1],[psMatPer{ii}(1:nl);0;zeros(nlmax-nl,1)],'LineWidth',2)
    hold on
%     plot(dr*[0:nl-1,nl-.5,nl:nlmax-1],dr*ones(1,nlmax+1),'LineWidth',2)
    set(gca,'ylim',[zmin,zmax],'xlim',[0 2],'FontSize',16)
    xlabel('$r/R_o$','interpreter','LaTeX','FontSize',16)
    ylabel('$\frac{P_sR_o}{2\sigma}\ \ \ \ \ $','interpreter','LaTeX','FontSize',24,'rotation',0)
    title(['$tU_0/R_o =\ $',num2str(tvecplot(ii)*U0/Ro)],'interpreter','LaTeX','FontSize',16)
    grid on
    hold off
    pause
end
plot(tvecplot*V0/Ro,f/dropWeight,'k','LineWidth',2)
hold on
plot(tvecplot*V0/Ro,numl(1:length(tvecplot)),'color',[.5 .5 .5],'LineWidth',2)

Integral = f*dtvecplot'
% ref = dropmass*g*1/40

hold on
grid on

numlMatPer = (psMatPer>0);
numlPer = sum(numlMatPer);
dnuml = [0,numlPer(2:end)-numlPer(1:end-1)];

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
set(gca,'xlim',[0 18],'ylim',[0 70],'FontName','Times','FontSize',24)%,'Xtick',[-.5:.25:.5],'Ytick',[]);
% saveas(gcf,['forcesOm0',num2str(10*Omega),'0Gamma',num2str(100*Gamma),'.fig'],'fig')
% print(fi,'-depsc','-r300',['Forces',num2str(10*Omega),'0Gamma',num2str(100*Gamma),'.eps'])



% figure
% for ii = size(etaMatPer,2)/2:15:size(etaMatPer,2)
%     plot([fliplr(-1*r),r(2:end)],(ii-1)*delta+[flipud(etaMatPer(:,ii));etaMatPer(2:end,ii)],'k','LineWidth',2)
%     grid on
%     hold on
% %     plot(r,zplot(ii)+zs)
% %     axis equal
%     set(gca,'xlim',[-.5 .5],'ylim',[.08 .19]);
% end






