%clc
%clear
close all;
addpath(fullfile(pwd, "..", "simulation_code"));
p = fullfile(pwd, "..",  "D50Quant100\rho1000sigma7220nu98muair0\RhoS1000SigmaS7220\R0350mm\ImpDefCornerAng180U39\N=20tol=5.00e-05"); % uigetdir();
cd(p);

global errored
errored = ~isfile('z.mat');
try

    load('ProblemConditions.mat'); NN = N;
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
load_vars('pressure_amplitudes.mat');
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

subplots = 5;
subplotWidths = 0.9;
subplotHeight = 0.85;
idxs = floor(linspace(1, size(etaMatPer,2)*0.715, subplots));
base = 150;
figure('Position', [5, 5, base*subplots/subplotWidth, 3*base/subplotHeight+200]);
for jj = 1:subplots
    ii = idxs(jj);

    thetaplot = linspace(0, pi, 100);%-%-0:thetaVec(end)/400:thetaVec(end);
    %-%-xsTop = xsoftheta(thetaplot,A2New,A3New);
    %-%-zsTop = zsoftheta(thetaplot,A2New,A3New);
    position = [(jj-1)*subplotWidths/subplots + (1-subplotWidths)/2, ...
        2/3 * subplotHeight + (1-subplotHeight)/2+0.03 , subplotWidths/subplots, subplotHeight/3];
    subplot('Position', position);
    
    zsTop = zs_from_spherical(thetaplot, oscillation_amplitudes(:, ii));
    xsTop = r_from_spherical(thetaplot, oscillation_amplitudes(:, ii)); 
    plot([-xsTop(end:-1:2), xsTop],[zsTop(end:-1:2), zsTop]+z(ii),'k','Linewidth',2);
    hold on

%%
    plot([fliplr(-1*r),r(2:end)],(zbplot(ii)+[0;flipud(etaMatPer(:,ii));etaMatPer(2:end,ii);0]),...
        'color',[.4 .4 .4],'LineWidth',2)
    hold on
    
    grid on
    l = 3;
    set(gca,'xlim',[-l l],'ylim',[-(l-1) (l+1)], 'Xtick',[-2 0 2],'FontName','Times','FontSize',14);
    %axis equal
    %xlabel('   $x/R_o$   ','interpreter','Latex','FontName','Times','FontSize',18)
    if jj == 1
        %ylabel('   $\frac{z}{R_o}\ \ \ $   ','interpreter','Latex','FontName','Times',...
        %    'FontSize',24,'rotation',0)
        %set(gca, 'Ytick', 0);
    else
        
    end
    yticklabels("");
        a = gca;
    a.XRuler.TickLabelGapOffset = -4; 
    
    

    title(sprintf("$ t/t_\\sigma =\\ $ %3.2f", tvec(ii)),'FontSize',18,...
            'interpreter','latex','FontName','Times')    
    drawnow
    %currFrame = getframe(gcf);
    %writeVideo(vidObj,currFrame);
    hold off
    
    % Pressure field distribution
    f = zeta_generator(pressure_amplitudes(:, ii));
    pfield = f(thetaplot) - sum(pressure_amplitudes(:, ii));
    
    position = [(jj-1)*subplotWidths/subplots + (1-subplotWidths)/2, ...
        1/3 * subplotHeight + (1-subplotHeight)/2, subplotWidths/subplots, subplotHeight/3];
    subplot('Position', position);
    plot(thetaplot*180/pi, pfield);
    set(gca,'xlim',[0 180], 'ylim', [-0.5, 2], 'Xtick', [45, 90, 135], ...
        'Ytick', [0 1 2], 'FontName','Times','FontSize',14);
    grid on
    %axis equal
    %xlabel('   $x/R_o$   ','interpreter','Latex','FontName','Times','FontSize',18)
    if jj == 1
        ylabel('$p(\theta)$ vs $\theta$','interpreter','Latex','FontName','Times',...
            'FontSize',24,'rotation',90)
    else
        yticklabels("");
        %set(gca,'ytick',[]);
    end
    
    a = gca;
    a.XRuler.TickLabelGapOffset = -4;

    drawnow
    hold off
    
    % Pressure amplitudes
    position = [(jj-1)*subplotWidths/subplots + (1-subplotWidths)/2, ...
        (1-subplotHeight)/2-0.04, subplotWidths/subplots, subplotHeight/3];
    subplot('Position', position);
    bar(0:NN, [-sum(pressure_amplitudes(:, ii)); pressure_amplitudes(:, ii)]);
    set(gca, 'ylim', [-0.7, 0.7], 'Xtick', [0, floor(NN/2)], 'Ytick', [-.5 0 .5], ...
        'FontName','Times','FontSize',14);
    grid on
    %axis equal
    %xlabel('   $x/R_o$   ','interpreter','Latex','FontName','Times','FontSize',18)
    if jj == 1
        ylabel('$B_{\ell}$ vs $\ell $','interpreter','Latex','FontName','Times',...
            'FontSize',24,'rotation',90)
    else
        yticklabels("");
        %set(gca,'ytick',[]);
    end
    a = gca;
    a.XRuler.TickLabelGapOffset = -4;
    drawnow
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

