%clc
close all;
%cd(fullfile(mfilename('fullpath'),  "..", ".."));
addpath(fullfile(pwd, "..", "simulation_code"));
%p = uigetdir();
curr = pwd;
saving_figure = fullfile(pwd, "..", "D50Quant100", "rho1000sigma7220nu98muair0", "RhoS1000SigmaS7220", "R0350mm", "ImpDefCornerAng180U39", "N=20tol=5.00e-05"); % uigetdir();
cd(saving_figure);

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
load_vars('numl.mat');
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

cd(saving_figure);

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
% set(vidObj,'FrameRate',10)
% open(vidObj);

%Gamma = 0; %%
%wzero = 1; %%
%thetaZero = 0; %%
%zb = Gamma/(Fr*wzero^2)*cos(wzero*tvec+thetaZero); %Elevation of the pool
%zbplot=zb; %(1:end-1);

saving_figure = figure();
saving_figure.Position = [100, 300, 866*0.75, 428*0.75];
N = floor(size(etaMatPer, 2)*0.8); M = floor(1.2/dr);
pfield_radial = zeros(M, N);
indexes = floor(linspace(1, size(etaMatPer,2), N));
plot_oscillations = oscillation_amplitudes(:, indexes);

hold on;

detatch_points = find(numl(1:end-1) ~= 0 & numl(2:end) == 0);
connect_points = find(numl(1:end-1) == 0 & numl(2:end) ~= 0);
%xline(tvec(detatch_points(1)), 'k--', 'LineWidth', 3, 'DisplayName', 'Contact ends');
for ii = 1:length(connect_points)
    dpname = '';
    if ii  == 1
        dpname = 'Contact region';
    end
    
    p = patch([tvec(connect_points(ii)) tvec(connect_points(ii)) tvec(detatch_points(ii)) tvec(detatch_points(ii))], ...
        [-.35 .35 .35 -.35], [234 182 118]/256, 'FaceAlpha', 0.15, 'EdgeColor', [.5 .5 .5], ...
        'DisplayName', dpname);
    if ii >= 1
        set(get(get(p, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    end
end

idxs = 2:20;
idxs = idxs(max(abs(plot_oscillations(1:idxs(end), :)), [] , 2) > 5e-3); %2.^(1:floor(log2(nb_harmonics)));
%cmap = colormap('spring'); %disp(size(cmap));
%numColors = size(cmap, 1);
lol = jet(length(idxs)+2);
colororder(flipud(lol(3:end, :)));%cmap(floor(linspace(1, numColors, length(idxs))), :));
%disp(size(times)); disp(size(deformation_amplitudes));
%plot(tvec(indexes), (plot_oscillations(idxs, :)), 'LineWidth',2);

for ii = 2:11
     plot(tvec(indexes), plot_oscillations(ii, :), 'LineWidth', 5-3*log10(ii),  'DisplayName', sprintf("A_{%d}(t)", ii));
end


% Create a figure
% Add labels and title
set(gca,'FontName','Times','FontSize',18);
set(gca,'defaultTextInterpreter','tex', ...
          'defaultAxesTickLabelInterpreter','tex', ...
          'defaultLegendInterpreter','tex');
xlabel('  t/T_d   ','FontName','Latin Modern Roman','FontSize',22, 'FontAngle', 'italic')
ylabel(' r/R_d','FontName','Latin Modern Roman',...
    'FontSize',22,'rotation',90, 'FontAngle', 'italic')
%legend(arrayfun(@(i) sprintf("Mode %d", i), idxs), 'FontSize', 15);
hl = legend('show', 'Location','eastoutside', 'FontSize', 14, 'FontName', 'Times');
yl = get(gca, 'YLim'); yl = [-.8*max(abs(yl)), .8*max(abs(yl))];
set(gca, 'YLim', yl); set(gca, 'XLim', [tvec(1) tvec(end)]);
box on
grid on
%xlim([0, 1]); ylim([0, 5]);
%title('Contact radius and pressure field evolution');

% Adjust axes
% === EXTRA LEGEND VIA GHOST AXIS (ON TOP) ================================
ax = gca;

% Ghost axis exactly over the main axis, invisible & non-interactive
axGhost = axes('Position', ax.Position, 'Color','none', 'Visible','off', ...
    'XLim', ax.XLim, 'YLim', ax.YLim, 'HitTest','off', ...
    'XColor','none','YColor','none');
if isprop(axGhost,'PickableParts'), set(axGhost,'PickableParts','none'); end

% Dummy swatch for "Contact region" (match your patch style)
hContact = patch('XData',NaN,'YData',NaN, ...  % NaNs so it doesn't affect limits
    'FaceColor',[234 182 118]/256,'FaceAlpha',0.15, ...
    'EdgeColor',[.5 .5 .5],'LineWidth',1, 'Parent',axGhost, ...
    'DisplayName','Contact region');

% (Optional) another dummy entry, e.g., detachment line style
% hDet = line(NaN,NaN,'LineStyle','--','Color','k','LineWidth',1.5, ...
%     'Parent',axGhost,'DisplayName','Detachment');

% Build the extra legend
lgd2 = legend(axGhost, hContact, {'Contact region'}, ...
    'Location','northwest', 'FontSize',18, 'FontName','Times', ...
    'Interpreter','tex', 'AutoUpdate','off');

% Make sure the legend is actually on top of everything
uistack(lgd2,'top');              % <-- this is the key fix
set(lgd2, 'Color','w');   % ghost legend
linkprop([ax axGhost], {'Position','XLim','YLim'});  % keep aligned
% ========================================================================

% ========================================================================

cd(curr);
set(gcf,'Renderer','painters');
saveas(saving_figure, "../../0_data/manual/amplitude_plotter_paper", 'fig');
%print(saving_figure, '-depsc', '-r300', "../../0_data/manual/amplitude_plotter_paper.eps");
exportgraphics(saving_figure, "../../0_data/manual/amplitude_plotter_paper.eps", 'ContentType','vector');


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



