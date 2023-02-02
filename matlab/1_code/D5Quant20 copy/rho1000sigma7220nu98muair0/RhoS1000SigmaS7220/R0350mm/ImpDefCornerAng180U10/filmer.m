clear
clc
close all

load('r.mat','r')
load('tExpLee.mat')
load('zExpLee.mat')
load('wzero.mat','wzero')
load('Ro.mat','Ro')
load('tvec.mat','tvec')
load('z.mat','z')
load('xdrop.mat','xdrop')
load('zdrop.mat','zdrop')
% load('phizero.mat','phizero')
phizero = 0;
load('Gamma.mat','Gamma')
load('g.mat','g')
% load('drop.mat','drop')
load('Ro.mat','Ro')
load('nr.mat','nr')
load('vz.mat')

Vo = abs(vz(1));
lfarad = 2*pi/12.6447;

tExpLee = tExpLee+.71302/1000;

zexp=zeros(size(tvec));
for j=1:length(tvec)
    if tvec(j) >= tExpLee(end)
        zexp(j) = 10;
    else
        indic = find(tExpLee>tvec(j),1)
        if indic == 1
            zexp(j) = z(j);
        else      
            zexp(j) = -Ro+(tExpLee(indic)-tvec(j))/(tExpLee(indic)-tExpLee(indic-1))*zExpLee(indic-1)+(tvec(j)-tExpLee(indic-1))/(tExpLee(indic)-tExpLee(indic-1))*zExpLee(indic);
        end
    end
%     pause
end


fig=figure

xplot=[-fliplr(r(2:nr+1)),r];

width = 300;

vidObj = VideoWriter(['FilmRo',num2str(Ro),'.mp4'],'MPEG-4')
set(vidObj,'FrameRate',24,'Quality',100);
open(vidObj);


for j=1:size(tvec,2)
    zp = -Gamma*g/(wzero^2)*sin(wzero*tvec(j)+phizero);
    load(['eta',num2str(j),'.mat'],'eta')
    etaplot = [flipud(eta(2:nr+1));eta];
    plot(xdrop/Ro,(zp+zexp(j)+zdrop)/Ro,'color',[.5 .5 .5],'Linewidth',4)
    hold on
    plot(xdrop/Ro,(zp+z(j)+zdrop)/Ro,'k','Linewidth',4)
    hold on
    plot(xplot(nr+1-width:nr+1+width)/Ro,(zp+etaplot(nr+1-width:nr+1+width))/Ro,'b','LineWidth',4);
    hold off
    axis equal
    axis([-5 5 -3.5 1.5])
    grid on
    xlabel('   $x/R_o$   ','interpreter','latex','FontSize',24,'interpreter','latex','FontName','Times')
    ylabel('   $\frac{z}{R_o}\ \ \ \ $      ','interpreter','latex','FontSize',36,'Rotation',0,'interpreter','latex','FontName','Times')
    set(gca,'FontSize',16,'XTick',[-5:5],'YTick',[-4:1],'FontName','Times')
%     if round(100*tvec(j)*80)/100 == round(tvec(j)*80)
%         numero = [num2str(round(100*tvec(j)*80)/100),'00'];
%     elseif round(100*tvec(j)*80)/100 == round(10*tvec(j)*80)/10
%         numero = [num2str(round(100*tvec(j)*80)/100),'0'];
%     else
%         numero = num2str(round(100*tvec(j)*80)/100);
%     end
    title(['$   tV_0/R_o = $',num2str(tvec(j)*Vo/Ro)],'FontSize',24,'interpreter','latex','FontName','Times')
%     title(['   t = ',num2str(1000*tvec(j)),' $[ms]$   '],'FontSize',24,'interpreter','latex')
    drawnow
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
%     if j==84
%         title('a','FontSize',24,'interpreter','latex','FontName','Times')
%         print('-depsc','-r300','a.eps')
%     elseif j==141
%         title('b','FontSize',24,'interpreter','latex','FontName','Times')
%         print('-depsc','-r300','b.eps')
%     elseif j==257
%         title('c','FontSize',24,'interpreter','latex','FontName','Times')
%         print('-depsc','-r300','c.eps')
%     elseif j==602
%         title('d','FontSize',24,'interpreter','latex','FontName','Times')
%         print('-depsc','-r300','d.eps')
%     end
end
close(vidObj);

