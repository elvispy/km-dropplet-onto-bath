clear
clc 
close all

load('etaS30n1.mat')
load('zS30n1.mat')
load('r.mat')
load('zdrop.mat')
load('xdrop.mat')
load('Ro.mat')


a=figure
plot(xdrop/Ro,(zdrop+zS30(1))/Ro,'k','LineWidth',1)
axis equal
hold on
plot(r(1:end-1)/Ro,etaS30(:,1)/Ro,'color',[.6 .6 .6],'LineWidth',1)
grid on
set(gca,'xlim',[0.5 .9],'ylim',[-.1 .2],'FontSize',24)
% xlabel(' $r/R_0$ ','interpreter','Latex','FontSize',20)
% ylabel(' $\frac{z}{R_0}\ $ ','interpreter','Latex','FontSize',24,'rotation',0)
title('   a   ','FontName','Times','FontSize',24)
print(a,'-depsc','-r300','a.eps')

b=figure
plot(xdrop/Ro,(zdrop+zS30(30))/Ro,'k','LineWidth',1)
axis equal
hold on
plot(r(1:end-1)/Ro,etaS30(:,30)/Ro,'color',[.6 .6 .6],'LineWidth',1)
grid on
set(gca,'xlim',[0.5 .9],'ylim',[-.1 .2],'FontSize',24)
% xlabel(' $r/R_0$ ','interpreter','Latex','FontSize',20)
% ylabel(' $\frac{z}{R_0}\ $ ','interpreter','Latex','FontSize',24,'rotation',0)
title('   b   ','FontName','Times','FontSize',24)
print(b,'-depsc','-r300','b.eps')

c=figure
plot(xdrop/Ro,(zdrop+zS30(31))/Ro,'k','LineWidth',1)
axis equal
hold on
plot(r(1:end-1)/Ro,etaS30(:,31)/Ro,'color',[.6 .6 .6],'LineWidth',1)
grid on
set(gca,'xlim',[0.5 .9],'ylim',[-.1 .2],'FontSize',24)
% xlabel(' $r/R_0$ ','interpreter','Latex','FontSize',20)
% ylabel(' $\frac{z}{R_0}\ $ ','interpreter','Latex','FontSize',24,'rotation',0)
title('   c   ','FontName','Times','FontSize',24)
print(c,'-depsc','-r300','c.eps')

