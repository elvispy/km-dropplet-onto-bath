close all;
c = pwd;
%p = uigetdir();
p = fullfile(pwd, 'D50Quant100/rho1000sigma7220nu98muair0/RhoS1000SigmaS7220/R0350mm/ImpDefCornerAng180U38');
cd(p);

try
    load('U0.mat');
    load('Ang.mat'); Cang = Ang * pi / 180;
    load('Fr.mat');
catch
    load('ProblemConditions.mat');
end

load('vz.mat'); Vo = abs(vz(1));
load('numl.mat','numl');
% files = dir(fullfile(pwd, "etaMatPer*.mat"));
% N = length(files);
% etaAux = [];
% for i = 1:N
%     load(files(i).name);
%     etaAux = [etaAux, etaMatPer];
% end

%zeta_generator2 = @zeta_generator;
%zeta_generator = zeta_generator2;
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
load('Ro.mat','Ro');

cd ..
load('rhoS.mat')
load('Ma.mat')
cd ..
cd ..
load('dr.mat')
load('IntMat.mat')

cd(p);


files = dir(fullfile(pwd, "psMat*.mat"));
N = length(files);
psAux = [];
for i = 1:N
    load(files(i).name);
    psAux = [psAux, psMatPer];
end
psMatPer = psAux; pss = psAux; save('ps.mat', 'pss');

f1 = figure(1);

hold on

% North Pole, South pole and bath
deep_blue = [66 145 245]/255;
verdinho = [0, 1, 0];
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
cd ..
cd 0_data/manual

FreeSurf = plot(tvec_p,etaOri(index_to_plot),'color',deep_blue,'LineWidth',4);

South = plot(tvec_p,z(index_to_plot)+south(index_to_plot),'color',verdinho,'LineWidth',2);
North = plot(tvec_p,z(index_to_plot)+north(index_to_plot),'color',verdinho,'LineWidth',2);

Center = plot(tvec_p,z(index_to_plot),'k','LineWidth',4);
set(gca,'FontSize',16); %,'xlim',[0 16],'ylim',[-2 8])
xlabel('   $t/t_s $   ','interpreter','LaTeX','FontSize',26)
ylabel('   $z/R \ \ \ \ $    ','interpreter','LaTeX','FontSize',26,'Rotation',0)
grid on
%ylim([0, 1]);
xlim([0, tvec_p(end)]);
saveas(f1,sprintf('CenterLineRadiusU0%g%.2fmm.fig',U0, 10*Ro),'fig')
print(f1,'-depsc','-r300',sprintf('CenterLineRadiusU0%g%.2fmm.eps',U0, 10*Ro))

% Pressed radius
f2 = figure(2);
PressedRad = plot(tvec_p,dr*numl(index_to_plot), 'color', deep_blue, 'LineWidth', 4);

set(gca,'FontSize',16); %,'xlim',[0 16],'ylim',[-2 8])
xlabel('   $t/t_s $   ','interpreter','LaTeX','FontSize',26)
ylabel('$r_c/R \ \ \ \ \ \ $    ','interpreter','LaTeX','FontSize',26,'Rotation',0)
ylim([0, 1]);
xlim([0, tvec_p(end)]);
grid on
saveas(f2,sprintf('PressedRadiusU0%g%.2fmm.fig',U0, 10*Ro),'fig')
print(f2,'-depsc','-r300',sprintf('PressedRadius%g%.2fmm.eps',U0, 10*Ro))

f3 = figure(3);

max_contact_radius = plot(tvec_p, max_width(index_to_plot), 'color', deep_blue, 'LineWidth', 4);
set(gca, 'FontSize', 16);
xlabel('  $  t/t_s $  ', 'interpreter', 'LaTeX', 'FontSize', 26);
ylabel(' .  $ w/R $ \ \ \  ', 'interpreter', 'LaTeX', 'FontSize', 26, 'Rotation', 0);
ylim([0.8, 1.2]);
xlim([0, tvec_p(end)]);
grid on
saveas(f3, sprintf("MaximumWidthU0%g%.2fmm.fig", U0, 10*Ro), 'fig');
print(f3, '-depsc', '-r300', sprintf('MaximumWidthU0%g%.2ffmm.eps', U0, 10*Ro));
%Integral = f*dtvecplot';
% daspect([.83 .73/5 1])

% %
% index1 = find(z<0,1);
% index2 = find(z(index1:end)>0,1);
% 
% tImpact = (tvec(index1)+tvec(index1-1))/2; save('tImpact.mat','tImpact')
% Uo = (vz(index1)+vz(index1-1))/2;save('Uo.mat','Uo')
% tend = (tvec(index1+index2-1)+tvec(index1+index2-2))/2;save('tend.mat','tend')
% Uend = (vz(index1+index2-1)+vz(index1+index2-2))/2;save('Uend.mat','Uend')
% tcont = tend-tImpact;save('tcont.mat','tcont')
% CRref = -Uend/Uo; save('CRref.mat','CRref')