%% This script plots 

%clc
clear vidObj
close all;
cd ..
addpath(fullfile(pwd, "simulation_code" ));
curr = pwd;
%p = uigetdir();
p = fullfile(pwd,  "D50Quant100", "rho1000sigma7220nu98muair0", "RhoS1000SigmaS7220", "R0350mm", "ImpDefCornerAng180U39", "N=20tol=5.00e-05");
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
load_vars("numl.mat");
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


cd ..
load('rho.mat','rho')

load('nu.mat','nu')
load('muair.mat')
load('g.mat','g') %gravitational constant

cd ..

load('nr.mat','nr')
load('dr.mat','dr')
load('r.mat')
load('zs.mat','zs')


%xplot = dr*(0:nr-1); save('xplot.mat','xplot')%I might remove or relocate this
load('xplot.mat')

cd(p);

ntimes = size(etaMatPer,2);

zplot = z(1:end-1);
etaOriplot = etaOri(1:end-1);

zmin = min(min(etaMatPer));
zmax = max(max(etaMatPer));

if errored
    file_name = 'errored_WavesAndDrop.mp4';
else
    file_name = 'WavesAndDrop.mp4';
end
%vidObj = VideoWriter(file_name,'MPEG-4');
%set(vidObj,'FrameRate',10)
%open(vidObj);

Gamma = 0; %%
wzero = 1; %%
thetaZero = 0; %%
zb = Gamma/(Fr*wzero^2)*cos(wzero*tvec+thetaZero); %Elevation of the pool
zbplot=zb; %(1:end-1);

saving_figure = figure();
saving_figure.Position = [100 100 800 400];

P_unit_old = M_unit / L_unit / T_unit^2;
P_unit = rhoS * U0^2;
thetaplot = linspace(0, pi, 200); jet2 = jet;
for ii = [60, 430, 715]%floor(linspace(1, size(etaMatPer,2), 300))
    %if tvec(ii)
    %% Drop video
%     subplot(1, 3, 2);
%     thetaplot = linspace(0, pi, 100);
%     zsTop = zs_from_spherical(thetaplot, oscillation_amplitudes(:, ii));
%     xsTop = r_from_spherical(thetaplot, oscillation_amplitudes(:, ii)); 
%     plot([-xsTop(end:-1:2), xsTop],[zsTop(end:-1:2), zsTop]+z(ii),'k','Linewidth',2);
%     hold on
% 
%     plot([fliplr(-1*r),r(2:end)],(zbplot(ii)+[0;flipud(etaMatPer(:,ii));etaMatPer(2:end,ii);0]),...
%         'color',[.4 .4 .4],'LineWidth',2)
%     hold on
%     %axis equal
%     grid on
%     set(gca,'xlim',[-2.5 2.5],'ylim',[-2 3],'Xtick',-5:5,'FontName','Times','FontSize',14);
%     xlabel('   $x/R_o$   ','interpreter','Latex','FontName','Times','FontSize',14)
%     ylabel('$z/R_o$','interpreter','Latex','FontName','Times',...
%         'FontSize',20,'rotation',90)
% %     t = (ii-1)/360;
% %     to = floor(t);
% %     t1 = floor(10*(t-to));
% %     t2 = round(100*(t-to-t1/10));
% %     if t2 == 10
% %         t2=0;
% %         t1=t1+1;
% %     end
% %     if t1 == 10
% %         t1 = 0;
% %         to = to+1;
% %     end
%     hold off
%     title(sprintf("$ t/t_\\sigma =\\ $ %3.2f", tvec(ii)),'FontSize',18,...
%             'interpreter','latex','FontName','Times')   

    %% Pressure field distribution
    hold on;
    f = zeta_generator(pressure_amplitudes(:, ii));
    pfield = f(thetaplot) - sum(pressure_amplitudes(:, ii));
    pmean = mean(pfield(1:50));
    pefield = P_unit_old / P_unit * pfield;
    pmean   = P_unit_old / P_unit * pmean;
    %position = [(jj-1)*subplotWidths/subplots + (1-subplotWidths)/2, ...
    %    1/3 * subplotHeight + (1-subplotHeight)/2, subplotWidths/subplots, subplotHeight/3];
    %subplot('Position', position);
    %subplot(1, 3, 1);
    plot(thetaplot*180/pi, pfield-pmean, 'DisplayName',sprintf("t = %.1f", tvec(ii)), ...
        'Color',jet2(mod(ii, 257), :), 'LineWidth',2);
    th = theta_from_cylindrical(numl(ii) * dr, oscillation_amplitudes(:, ii));
    xline(th*180/pi, '--', 'HandleVisibility','off', 'Color',jet2(mod(ii, 257), :), ...
        'LineWidth',2);
    
end
set(gca,'xlim',[0 180], 'ylim', [-0.5, 2], 'Xtick', [0, 45, 90, 135, 180], ...
        'Ytick', [0 1 2], 'FontName','Times','FontSize',14);
grid on
legend show
legend('Location','northwest', 'FontSize',18);
%axis equal
%xlabel('   $x/R_o$   ','interpreter','Latex','FontName','Times','FontSize',18)
xlabel('$\theta$','interpreter','Latex','FontName','Times','FontSize',14); 
ylabel('$p(\theta)/(\rho V_o^2)$','interpreter','Latex','FontName','Times',...
    'FontSize',20,'rotation',90)


a = gca;
a.XRuler.TickLabelGapOffset = -4;

cd(curr);
saveas(saving_figure, "../0_data/manual/pressure_dist", 'fig');
print(saving_figure, '-depsc', '-r300', "../0_data/manual/pressure_dist.eps");


    %% Pressure amplitudes
    
    %position = [(jj-1)*subplotWidths/subplots + (1-subplotWidths)/2, ...
    %    (1-subplotHeight)/2-0.04, subplotWidths/subplots, subplotHeight/3];
    % subplot(1, 3, 3);
    % bar(0:NN, [-sum(pressure_amplitudes(:, ii)); pressure_amplitudes(:, ii)]);
    % set(gca, 'ylim', [-0.7, 0.7], 'Xtick', [0, floor(NN/2)], 'Ytick', [-.5 0 .5], ...
    %     'FontName','Times','FontSize',14);
    % grid on
    % xlabel('$\ell$','interpreter','Latex','FontName','Times','FontSize',14); 
    % ylabel('$B_{\ell}$','interpreter','Latex','FontName','Times',...
    %     'FontSize',20,'rotation',90)
    % 
    % a = gca;
    % a.XRuler.TickLabelGapOffset = -4;
    % 
    % drawnow
    % %g = gcf;
    % %g.WindowState = 'maximized';
    % currFrame = getframe(h);
    % writeVideo(vidObj,currFrame);
    % hold off
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



