clc
clear
close all

load('r.mat')
load('zs.mat','zs')
load('Gamma.mat')
load('wzero.mat')
load('g.mat')
load('thetaZero.mat')
load('xdrop.mat')
load('zdrop.mat')
load('Ro.mat')

load('z.mat')
load('etaOri.mat')
load('tvec.mat')

load('tExpLee.mat')
load('zExpLee.mat')

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

ntimes = size(etaMatPer,2);

zplot = z;
etaOriplot = etaOri;

zmin = min(min(etaMatPer));
zmax = max(max(etaMatPer));
caxis([zmin zmax]);

vidObj = VideoWriter('WavesAndDrop3D21paperExp.mp4','MPEG-4');
set(vidObj,'Quality',100,'FrameRate',100)
open(vidObj);

zb = g*Gamma/(wzero^2)*cos(wzero*tvec+thetaZero);
zbplot=zb;

thetaplot = 0:pi/180:2*pi;
[rsurf,thetasurf] = meshgrid(r,thetaplot);

thetadrop = 0:pi/180:pi;
[rdropsurf,thetadropsurf] = meshgrid(xdrop,thetadrop);
[zdropsurf,~] = meshgrid(zdrop,thetadrop);

for ii = 1:1:size(etaMatPer,2)
    [zsurf,~] = meshgrid(etaMatPer(:,ii),thetaplot);
    u = surf(rdropsurf.*cos(thetadropsurf),rdropsurf.*sin(thetadropsurf),zbplot(ii)+zexp(ii)+zdropsurf,...
        'FaceColor',[.5 .5 .5],'LineStyle','none')
    hold on
    w = surf(rsurf.*cos(thetasurf),rsurf.*sin(thetasurf),zbplot(ii)+zsurf,...
        'FaceColor',[.15 .5 .95],'LineStyle','none')%[.5 .5 1]
%     contour(rsurf.*cos(thetasurf),rsurf.*sin(thetasurf),zbplot(ii)+zsurf)
    axis equal
    axis off
   
    d = surf(rdropsurf.*cos(thetadropsurf),rdropsurf.*sin(thetadropsurf),zbplot(ii)+zplot(ii)+zdropsurf,...
        'LineStyle','none','FaceColor',[1 1 1])
%     alpha(d,.3)
%     hold on
%     grid on
%     set(gca,'xlim',[-1 1],'ylim',[-.2 .5],'FontName','Times','FontSize',24);
%     xlabel('   x[mm]   ','FontName','Times','FontSize',24)
%     ylabel('   y [mm]   ','FontName','Times','FontSize',24)
    t = (ii-1)/360;
    to = floor(t);
    t1 = floor(10*(t-to));
    t2 = round(100*(t-to-t1/10));
    if t2 == 10
        t2=0;
        t1=t1+1;
    end
    if t1 == 10
        t1 = 0;
        to = to+1;
    end
    title(['   t/T_f = ',num2str(to),'.',num2str(t1),num2str(t2)],'FontName','Times','FontSize',24)
    
    luz  = light('Position',10*[0   12 .2],'style','local')
%     luz1 = light('Position',[-10   0 .1])%,'style','local')
    luz2 = light('Position',[-10 0 10])%,'style','local')
%     luz3 = light('Position',.2*[-12 -.3 .4],'style','local')
%     luz4 = light('Position',.2*[-12 -.4 .4],'style','local')
%     luz5 = light('Position',.2*[-12 -.5 .4],'style','local')
    material([.2 0 1])
%     Camera([5, 5, 5], [0, 0, 0], PI/4):
    campos(4*[10 0 7])
    camtarget([0 0 .05])
    camzoom(1)
    camva(pi)
%   7
    
    
    
    drawnow
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
    hold off

end
close(vidObj);






