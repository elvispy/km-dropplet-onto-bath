close all;
c = pwd;
disp(mfilename('fullpath'));
%p = uigetdir();
addpath(fullfile(pwd, "..", "simulation_code" ));
p = fullfile(pwd, "..",  "D50Quant100", "rho1000sigma7220nu98muair0", "RhoS1000SigmaS7220", "R0350mm", "ImpDefCornerAng180U39", "N=20tol=5.00e-05");
%p = uigetdir();
cd(p);

savefile = true; % flag to export files

try
    load('ProblemConditions.mat'); %Cang = Ang * pi / 180;
catch
    load('U0.mat');
    load('Fr.mat');
    load('Ang.mat'); %Cang = Ang * pi / 180;
    disp("Couldn't find Problem Conditions");
end

load('vz.mat'); Vo = abs(vz(1));
load('numl.mat','numl');

load('etas.mat');
load('z.mat');
load('etaOri.mat')
load('tvec.mat')
%load('Fr.mat')
load('oscillation_amplitudes.mat');
load('Rv.mat')
%load('dtb.mat');
% Rv = zeros(1, size(oscillation_amplitudes, 2));
% for ii = 1:size(oscillation_amplitudes, 2)
%     Rv(ii) = zs_from_spherical(pi, oscillation_amplitudes(:, ii));
% end
cd ..
cd ..
load('Ro.mat','Ro');

cd ..
load('rhoS.mat')
load('Ma.mat')
cd ..
cd ..
load('dr.mat')
load('IntMat.mat')

cd(p);


%files = dir(fullfile(pwd, "psMat*.mat"));
%N = length(files);
%psAux = [];
%for i = 1:N
%    load(files(i).name);
%    psAux = [psAux, psMatPer];
%end
%psMatPer = psAux; pss = psAux; save('ps.mat', 'pss');

%% Free surface and center
f1 = figure(1);
f1.Position = [50 457 560*0.75 420*0.75];
hold on

deep_blue = [13 120 245]/255;
verdinho = [0, .6, 0];
index_to_plot = 1:floor(length(z)*0.8);
tvec_p = tvec(index_to_plot);
south = zeros(1, size(oscillation_amplitudes, 2));
north = zeros(1, size(oscillation_amplitudes, 2));
for ii = 1:size(oscillation_amplitudes, 2)
    south(ii) = zs_from_spherical(pi, oscillation_amplitudes(:, ii));
    north(ii) = zs_from_spherical(0, oscillation_amplitudes(:, ii));
end

max_width = zeros(1, size(oscillation_amplitudes, 2));
for ii = 1:size(oscillation_amplitudes, 2)
    max_width(ii) = maximum_contact_radius(oscillation_amplitudes(:, ii));
end

cd(c);
%cd ..
cd("../../0_data/manual");

FreeSurf = plot(tvec_p,etaOri(index_to_plot),'color',deep_blue,'LineWidth',4, 'LineStyle', '--');

hold on;  plotter('surface'); plotter('center_');
Center = plot(tvec_p,z(index_to_plot),'k--','LineWidth',4);
set(gca,'FontSize',16); %,'xlim',[0 16],'ylim',[-2 8])
xlabel('   $t/T_s $   ','interpreter','LaTeX','FontSize',26)
ylabel('$z/R_s$','interpreter','LaTeX','FontSize',26,'Rotation',90)
text(-1,2.1,"(a)", 'FontSize', 20);
grid on
%ylim([0, 1]);
xlim([0, tvec_p(end)]);


%% South and north pole
%f4 = figure(4);
%f4.Position = [13 14 560*0.75 420*0.75];
%hold on;
South = plot(tvec_p,z(index_to_plot)+south(index_to_plot),'color',verdinho,'LineWidth',2);
North = plot(tvec_p,z(index_to_plot)+north(index_to_plot),'color',verdinho,'LineWidth',2);
plotter('bottom_'); plotter('top_');


if savefile
    saveas(f1,sprintf('CenterLineRadiusU0%g%.2fmm.fig',U0, 10*Ro),'fig')
    print(f1,'-depsc','-r300',sprintf('CenterLineRadiusU0%g%.2fmm.eps',U0, 10*Ro))
end

%% Contact radius
f2 = figure(2);
f2.Position = [526 457 560*0.75 420*0.75];
PressedRad = plot(tvec_p,dr*numl(index_to_plot), 'color', deep_blue, 'LineWidth', 4);
hold on;  plotter('c_radius');
set(gca,'FontSize',16); %,'xlim',[0 16],'ylim',[-2 8])
xlabel('   $t/T_s $   ','interpreter','LaTeX','FontSize',26)
ylabel('$r_c/R_s$','interpreter','LaTeX','FontSize',26,'Rotation',90)
ylim([0, 1]);
xlim([0, tvec_p(end)]);
text(-1,.9,"(b)", 'FontSize', 20);
grid on
if savefile
    saveas(f2,sprintf('PressedRadiusU0%g%.2fmm.fig',U0, 10*Ro),'fig')
    print(f2,'-depsc','-r300',sprintf('PressedRadius%g%.2fmm.eps',U0, 10*Ro))
end

%% Max width evolution
f3 = figure(3);
f3.Position = [1006 457 560*0.75 420*0.75];
max_contact_radius = plot(tvec_p, max_width(index_to_plot), 'color', deep_blue, 'LineWidth', 4);
hold on;  plotter('width');
set(gca, 'FontSize', 16);
xlabel('  $  t/T_s $  ', 'interpreter', 'LaTeX', 'FontSize', 26);
ylabel('$ w/R_s $', 'interpreter', 'LaTeX', 'FontSize', 26, 'Rotation', 90);
ylim([0.8, 1.2]);
xlim([0, tvec_p(end)]);
text(-1,1.15,"(c)", 'FontSize', 20);
grid on
if savefile
    saveas(f3, sprintf("MaximumWidthU0%g%.2fmm.fig", U0, 10*Ro), 'fig');
    print(f3, '-depsc', '-r300', sprintf('MaximumWidthU0%g%.2ffmm.eps', U0, 10*Ro));
end
%Integral = f*dtvecplot';
% daspect([.83 .73/5 1])

cd(c);

% Extract all files with one linewidth 
function plotter(prefixString)
    l = 0; r = 5.1;
    % Get a list of all CSV files in the directory
    files = dir(fullfile("../../**/*.csv"));
    
    % Loop through each file
    for i = 1:length(files)
        % Get the file name
        fileName = files(i).name;
        
        % Check if the file name contains the prefix string
        if startsWith(fileName, prefixString)
            % Read the CSV file (assuming no headers)
            data = sortrows(readmatrix(fullfile(files(i).folder, fileName)), 1);
            data = data(data(:, 1) > l & data(:, 1) <= r, :);
            % Plot the data
            %figure;
            if contains(fileName, 'DNS')
                plot(data(:, 1), data(:, 2), 'k-');  % Adjust columns as needed
            elseif contains(fileName, '1PKM')
                plot(data(:, 1), data(:, 2), 'b-');
            else
                plot(data(:, 1), data(:, 2), 'r', 'LineWidth', 3);
            end
            %title(['Plot of ', fileName]);
            %xlabel('X-axis');
            %ylabel('Y-axis');
        end
    end
end