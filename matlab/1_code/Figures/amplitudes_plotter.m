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

Gamma = 0; %%
wzero = 1; %%
thetaZero = 0; %%
zb = Gamma/(Fr*wzero^2)*cos(wzero*tvec+thetaZero); %Elevation of the pool
zbplot=zb; %(1:end-1);

saving_figure = figure();
saving_figure.Position = [100 100 600 400];
N = floor(size(etaMatPer, 2)*0.8); M = floor(1.2/dr);
pfield_radial = zeros(M, N);
indexes = floor(linspace(1, size(etaMatPer,2), N));
plot_oscillations = oscillation_amplitudes(:, indexes);
hold on;
for ii = 2:size(oscillation_amplitudes, 1)
    plot(tvec(indexes), plot_oscillations(ii, :), 'DisplayName',num2str(ii));
    %jj = indexes(ii);
    %f = zeta_generator(oscillation_amplitudes(:, jj));

    %thetaplots = theta_from_cylindrical(dr .* (0:numl(jj)), oscillation_amplitudes(:, jj));
    %pfield_now = f(thetaplots) - sum(pressure_amplitudes(:, jj));
    %pmean = mean(f(linspace(0, pi/10, 20)) - sum(pressure_amplitudes(:, jj)));
    
    %pfield_radial(1:(numl(jj)+1), ii) = pfield_now'-pmean;


end

% Create a figure

% Plot the matrix using pcolor
%NN = ceil(1/dr);
%pfield_radial(pfield_radial > 2) = 2; 
%p = pcolor(dr * (0:(M-1)), tvec(indexes) , pfield_radial');
%hold on
%plot(numl(indexes) * dr, tvec(indexes), 'r--', 'LineWidth',1)
% Enhance the plot
%p.EdgeColor = 'none';            % Remove gridlines for a smoother look
%jet2 = jet; jet2(1, :) = 1;
%colormap(flipud(bone));                   % Choose a visually striking colormap (jet, parula, etc.)
%cb = colorbar;                        % Add a colorbar for reference
%shading interp;                  % Interpolate shading for a smoother gradient
%ylabel(cb, '$p/p_0$','interpreter','Latex','FontName','Times',...
%    'FontSize',20,'rotation',90)
% Add labels and title
set(gca,'FontName','Times','FontSize',14);
xlabel('  $t/t_{\sigma}$   ','interpreter','Latex','FontName','Times','FontSize',20)
ylabel(' $x/R_o$','interpreter','Latex','FontName','Times',...
    'FontSize',20,'rotation',90)
legend()
%xlim([0, 1]); ylim([0, 5]);
%title('Contact radius and pressure field evolution');

% Adjust axes
%axis tight;                      % Fit the plot closely around the data

cd(curr);
saveas(saving_figure, "../../0_data/manual/amplitude_plotter", 'fig');
print(saving_figure, '-depsc', '-r300', "../../0_data/manual/amplitude_plotter.eps");


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



