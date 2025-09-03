% This script will plot a number of variables relating such as coefficient
% of restitution, contact time and max deflection.


clear data; close all;
data.D = 50;
data.Quant = 50;
%rho = 1; % must multiply by x1000

%data.nu = 9.78E-3; % Multiply by x10000
%muair = 0;
%RhoS = 1; % must multiply by x1000
%data.sigma = 72.20; % must multiply by x100
%data.Ro = 0.035; % linspace(0.02, 0.05, 5)'; % must multiply by x10
%Oil
%data.Bo = 0.056;
%data.Oh = 0.058;

% Water
%data.Bo = 0.017;
%data.Oh = 0.006;
%Ang = 180;
%U = 18; %linspace(28, 50, 5)';
current_folder = pwd;
cd ..
files = dir("**/simulation_postprocessing.mat");
cd(current_folder);
Westar = []; Bo = []; Oh = []; max_deflection = []; contact_time = []; 
coef_restitution = []; N = [];
plotting_data = table(Westar, Bo, Oh, max_deflection, contact_time, coef_restitution, N);
standalone = true;
saving = false;

for ii = 1:length(files)
    
    try 
        if length(files) == 1; folder_name = files.folder; else; folder_name = files(ii).folder; end
        %if isfield(data, 'D') && isfield(data, 'Quant') && ~contains(folder_name, sprintf("D%dQuant%d", data.D, data.Quant)); continue; end
        cd(folder_name);
        try
            load("ProblemConditions.mat");
        catch
            load("U0.mat");
            load("T.mat");
        end
        cd ..
        if ~isfile("RoMaker.m"); cd ..; end
        
        if ~isfile("Ro.mat") && isfile("RoMaker.m"); RoMaker; end
            
        load('Ro.mat','Ro'); 
        
        cd ..
        %load('rhoS.mat','rhoS')%Sphere density
        %load('sigmaS.mat')%Sphere's surface tension
        cd ..
        if ~isfile("rho.mat") && isfile("BathMaker.m"); BathMaker; end
        load('rho.mat','rho')
        load('sigma.mat','sigma')
        load('nu.mat','nu')
        load('muair.mat', 'muair')
        load('g.mat','g') %gravitational constant
        if nu <=0; continue; end
        cd(folder_name);
        Westar = rho * U0.^2 * Ro / sigma;
        Bo = rho * g * Ro.^2 / sigma;
        Oh = nu / sqrt(sigma * Ro * rho);
        t_sigma = sqrt(rho * Ro^3/sigma);

        simul.Ro = Ro; simul.rho = rho; simul.sigma = sigma; simul.U0 = U0;
        simul.nu = nu; simul.We = Westar; simul.Bo = Bo; simul.Oh = Oh;
        simul.folder = folder_name; simul.N = N;

        if is_valid(simul, data)
            load('oscillation_amplitudes.mat');
            if norm(oscillation_amplitudes) < 0.01
                continue
            end
            load("simulation_postprocessing.mat");     
            if ~exist('N', 'var'); N = nan; end
            max_deflection = abs(max_def); if isempty(max_deflection) == true; max_deflection = NaN; end
            coef_restitution = CRref; if isempty(coef_restitution) == true; coef_restitution = NaN; end
            contact_time = tcont * T /t_sigma; if isempty(contact_time) == true; contact_time = NaN; end
            plotting_data = [plotting_data; {Westar, Bo, Oh, max_deflection, contact_time, coef_restitution, N}];
        end
    catch ME
        warning(ME.message);
        disp(pwd)
    end

end

% Just for  the comparison against chase
plotting_data = plotting_data(plotting_data.N == 60, :);


% Load experimental data from the Excel file
filename = 'DataForReboundsOnly.xlsx'; % Replace with the correct file path
cd(current_folder);
dataExp = readtable(filename, 'Sheet', 'Sheet1', 'Range', 'A:T'); % Adjust the range if necessary

% Extract the experimental data (We, Bo, Oh)
We_exp = str2double(dataExp{3:end, 17}); % Column 'We' from row 3 onward (index 14 corresponds to column 13 'We')
Bo_exp = str2double(dataExp{3:end, 18}); % Column 'Bo' from row 3 onward (index 15 corresponds to column 14 'Bo')
Oh_exp = str2double(dataExp{3:end, 19}); % Column 'Oh' from row 3 onward (index 16 corresponds to column 15 'Oh')
coef_res_exp = str2double(dataExp{3:end, 16});

% Clean up data by removing NaN values
valid_exp_data = ~isnan(Bo_exp) & ~isnan(Oh_exp);
We_exp = We_exp(valid_exp_data);
Bo_exp = Bo_exp(valid_exp_data);
Oh_exp = Oh_exp(valid_exp_data);
coef_res_exp  = coef_res_exp(valid_exp_data);
% Assuming 'plotting_data' table contains your simulation data
% Extract simulation data from the 'plotting_data' table
Westar_sim = plotting_data.Westar;
coef_restitution_sim = plotting_data.coef_restitution;
max_deflection_sim = plotting_data.max_deflection;
Oh_sim = plotting_data.Oh;
Bo_sim = plotting_data.Bo;

% Create figures for coef_restitution and max_deflection

% 1. Plot for Coefficient of Restitution (Simulation + Experimental Data)
figure;
scatter(Westar_sim, coef_restitution_sim, 150, Bo_sim, 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Simulation');
hold on;
scatter(We_exp, coef_res_exp, 30, Bo_exp, '^', 'filled', 'DisplayName', 'Experimental');
% Customize the plot for coef_restitution
set(gca, 'FontSize', 16, 'XScale', 'log');
xlabel('   $We$   ', 'interpreter', 'LaTeX', 'FontSize', 20);
ylabel('   $\alpha \ \ \ $    ', 'interpreter', 'LaTeX', 'FontSize', 20, 'Rotation', 0);
title('Coefficient of Restitution vs. Westar for varying Bo and Oh', 'interpreter', 'LaTeX', 'FontSize', 20);
colorbar; % Show color bar for Oh
colormap jet; % Apply the color map based on Oh values
legend('show');
set(gcf, 'Position', [1288 557 560 420]); % Adjust the figure size if needed
hold off;

% 2. Plot for Maximum Deflection (Simulation + Experimental Data)
figure;
scatter(Westar_sim, max_deflection_sim, 100, Bo_sim, 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Simulation');
hold on;
%scatter(We_exp, Bo_exp, 100, Oh_exp, '^', 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Experimental');

% Customize the plot for max_deflection
set(gca, 'FontSize', 16, 'XScale', 'log');
xlabel('   $We$   ', 'interpreter', 'LaTeX', 'FontSize', 20);
ylabel('   $z / R_o \ \ \ $    ', 'interpreter', 'LaTeX', 'FontSize', 20, 'Rotation', 0);
title('Maximum Deflection vs. Westar for varying Bo and Oh', 'interpreter', 'LaTeX', 'FontSize', 20);
colorbar; % Show color bar for Oh
colormap jet; % Apply the color map based on Oh values
legend('show');
set(gcf, 'Position', [680 557 560 420]); % Adjust the figure size if needed
hold off;


% Prepare your data (assuming plotting_data is your table)
% Example: plotting_data = table(Westar, Bo, Oh, contact_time, coef_restitution, max_deflection, N);

% Create a list of unique Bo values for custom marker styles
unique_Bo = unique(plotting_data.Bo);
markers = {'o', 's', '^', 'd', 'p', 'h', 'x', '+'}; % You can add or remove marker styles

% Create a colormap based on Oh
cmap = jet(length(unique(plotting_data.Oh)));

% Prepare the figure for contact_time plot
figure;
hold on;
grid on;

% Loop through the unique Bo values and plot each with a different marker style
for i = 1:length(unique_Bo)
    % Select the rows with the current Bo value
    idx = plotting_data.Bo == unique_Bo(i);
    
    % Scatter plot for the current subset of data
    fprintf("Bo = %g, min We = %.2e, max We = %.2e \n", unique_Bo(i), ...
        min(plotting_data.Westar(idx)), max(plotting_data.Westar(idx)));
    scatter(plotting_data.Westar(idx), plotting_data.contact_time(idx), ...
        100, plotting_data.Oh(idx), markers{mod(i-1, length(markers)) + 1}, 'filled', 'DisplayName', sprintf('Bo = %.2g', unique_Bo(i)));
end

% Customize the plot
set(gca, 'FontSize', 16, 'XScale', 'log');
xlabel('   $We$   ', 'interpreter', 'LaTeX', 'FontSize', 20);
ylabel('   $t^* \ \ \ $   ', 'interpreter', 'LaTeX', 'FontSize', 20, 'Rotation', 0);
title('Contact Time vs. Westar for varying Bo and Oh', 'interpreter', 'LaTeX', 'FontSize', 20);
colorbar; % Show color bar for Oh
colormap(cmap); % Apply the color map based on Oh values
legend('show');
set(gcf, 'Position', [98 557 560 420]); % Adjust the figure size if needed
hold off;



function bool = is_valid(simul, data)
    bool = true; tolerance = 0.1;
    fnames = fieldnames(data);
    for ii = 1:length(fnames)
        fieldname = fnames{ii};
        if fieldname == "Quant"
            % checking the name of the folder
            
            check = false;
            for folder = data.(fieldname)
                check = (check || ~isempty(regexp(simul.folder, ...
                sprintf("Quant%g", folder), 'ONCE')));
            end
        elseif fieldname == "D"
            % checking the name of the folder
            
            check = false;
            for folder = data.(fieldname)
                check = (check || ~isempty(regexp(simul.folder, ...
                sprintf("D%gQuant", folder), 'ONCE')));
            end
        elseif length(data.(fieldname)) < 2
            % checking that data is roughly equal to value
            check = relerr(data.(fieldname), simul.(fieldname)) < tolerance;
            
        elseif length(data.(fieldname)) == 2
            % Checking that data is within bounds
            list = data.(fieldname);
            check = (list(1) * (1-tolerance) < simul.(fieldname)) && ...
                    simul.(fieldname) < list(2) * (1 + tolerance);
        else
            warning("There is some field name in the checking data that does not conform to these rules");
            check = true;
        end
        bool = bool && check;
    end
end
function res = relerr(a, b)
    res = norm(a-b, 1)/max(abs([a, b]));
end
