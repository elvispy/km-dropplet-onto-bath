%clc
%clear
close all;
addpath(fullfile(pwd, "..", "simulation_code"));
p = fullfile(pwd, "..",  "D50Quant100\rho1000sigma7220nu98muair0\RhoS1000SigmaS7220\R0350mm\ImpDefCornerAng180U39\N=20tol=5.00e-05"); % uigetdir();
cd(p);

global errored
errored = ~isfile('z.mat');
try

    load('ProblemConditions.mat');
catch
    load('U0.mat');
    load('Fr.mat');
    disp("Couldn't find Problem Conditions");
end
load_vars('vz.mat'); Vo = abs(vz(1));


try
    load_vars('etas.mat');
    etaMatPer = etas;
catch
    files = dir(fullfile(pwd, "etaMatPer*.mat"));
    N = length(files);
    etaAux = [];
    for i = 1:N
        load_vars(files(i).name);
        etaAux = [etaAux, etaMatPer];
    end
    etaMatPer = etaAux;
end
load_vars('z.mat')
load_vars('etaOri.mat')
load_vars('tvec.mat')

load_vars('oscillation_amplitudes.mat');
Rv = zeros(1, size(oscillation_amplitudes, 2));
for ii = 1:size(oscillation_amplitudes, 2)
    Rv(ii) = zs_from_spherical(pi, oscillation_amplitudes(:, ii));
end

%%

cd ..
try
    load('Ro.mat','Ro')%Sphere's radius in CGS
catch
    cd ..
    load('Ro.mat','Ro')
end
cd ..
load('rhoS.mat','rhoS')%Sphere density
%load('sigmaS.mat')%Sphere's surface tension

cd ..
load('rho.mat','rho')
%load('sigma.mat','sigma')
load('nu.mat','nu')
load('muair.mat')
load('g.mat','g') %gravitational constant

cd ..
%load('D.mat')%Domain diameter in units of droplet radii
%load('quant.mat')%number of dr's contained in an undeformed dropelt radius
load('nr.mat','nr')
load('dr.mat','dr')
load('r.mat')
load('zs.mat','zs')
%load('Deload('zs.mat','zs')lta.mat','Delta')
%load('IntMat.mat','IntMat')
%load(sprintf('DTNnew345nr%dD%drefp10.mat', nr, D),'DTNnew345')
%DTN = DTNnew345;
%clear DTNnew345

%xplot = dr*(0:nr-1); save('xplot.mat','xplot')%I might remove or relocate this
load('xplot.mat')

%cd(['rho',num2str(1000*rho),'sigma',num2str(round(100*sigma)),'nu',num2str(round(10000*nu)),'muair',num2str(muair)])

%cd(['RhoS',num2str(rhoS*1000),'SigmaS',num2str(round(100*sigmaS))])
%load('Ma.mat','Ma')%Dimensionless mass of sphere
%load('Ra.mat','Ra')%Density ratio

%%
% cd ..
% load('Ro.mat');
% 
% cd ..
% load('rhoS.mat');

%cd ..
%load('r.mat')
%load('nr.mat')
%load('dr.mat')
%load('zs.mat','zs')
%load('Gamma.mat')
%load('wzero.mat')
%load('thetaZero.mat')
%load('xplot.mat')

%cd(['Rho',num2str(rhoS*1000)])

%cd(['ImpDefAng',num2str(round((pi-Cang)*180/pi)),'U',num2str(U0)])
%load('Rv.mat')
%Rv = ones(size(etaMatPer, 2), 1);
%load('Fr.mat')



cd(p);


ntimes = size(etaMatPer,2);

zplot = z(1:end-1);
etaOriplot = etaOri(1:end-1);

zmin = min(min(etaMatPer));
zmax = max(max(etaMatPer));

% if errored
%     file_name = 'errored_WavesAndDrop.mp4';
% else
%     file_name = 'WavesAndDrop.mp4';
% end
% vidObj = VideoWriter(file_name,'MPEG-4');
% set(vidObj,'FrameRate',20)
% open(vidObj);

Gamma = 0; %%
wzero = 1; %%
thetaZero = 0; %%
zb = Gamma/(Fr*wzero^2)*cos(wzero*tvec+thetaZero); %Elevation of the pool
zbplot=zb; %(1:end-1);

subplots = 7;
figure('Position', [100, 100, 1024, 400]);
subplotWidths = 0.8;
subplotHeight = 0.8;
for ii = floor(linspace(1, size(etaMatPer,2), subplots))


    thetaplot = linspace(0, pi, 100);%-%-0:thetaVec(end)/400:thetaVec(end);
    %-%-xsTop = xsoftheta(thetaplot,A2New,A3New);
    %-%-zsTop = zsoftheta(thetaplot,A2New,A3New);
    position = [(i-1)*subplotWidths/subplots + (1-subplotWidths)/2, (1-subplotHeight)/2, subplotWidths/subplots, subplotHeight];
    subplot('Position', position);
    
    zsTop = zs_from_spherical(thetaplot, oscillation_amplitudes(:, ii));
    xsTop = r_from_spherical(thetaplot, oscillation_amplitudes(:, ii)); 
    plot([-xsTop(end:-1:2), xsTop],[zsTop(end:-1:2), zsTop]+z(ii),'k','Linewidth',2);
    hold on

%%
    plot([fliplr(-1*r),r(2:end)],(zbplot(ii)+[0;flipud(etaMatPer(:,ii));etaMatPer(2:end,ii);0]),...
        'color',[.4 .4 .4],'LineWidth',2)
    hold on
    axis equal
    grid on
    set(gca,'xlim',[-5 5],'ylim',[-1.5 3],'Xtick',-5:5,'FontName','Times','FontSize',18);
    xlabel('   $x/R_o$   ','interpreter','Latex','FontName','Times','FontSize',18)
    ylabel('   $\frac{z}{R_o}\ \ \ $   ','interpreter','Latex','FontName','Times',...
        'FontSize',24,'rotation',0)
%     t = (ii-1)/360;
%     to = floor(t);
%     t1 = floor(10*(t-to));
%     t2 = round(100*(t-to-t1/10));
%     if t2 == 10
%         t2=0;
%         t1=t1+1;
%     end
%     if t1 == 10
%         t1 = 0;
%         to = to+1;
%     end

    title(sprintf("$ t/t_\\sigma =\\ $ %3.2f", tvec(ii)),'FontSize',18,...
            'interpreter','latex','FontName','Times')    
    drawnow
    %currFrame = getframe(gcf);
    %writeVideo(vidObj,currFrame);
    hold off
end
%close(vidObj);


function load_vars(str)
    global errored
    
    if errored == true
        str = "errored_" + str; 
    end
    
    vars = load(str);
    fn = fieldnames(vars);
    for ii = 1:length(fn)
        assignin('caller', fn{ii}, vars.(fn{ii}));
    end
    
end



