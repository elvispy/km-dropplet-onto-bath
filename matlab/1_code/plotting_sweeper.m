% This script will plot a number of variables relating such as coefficient
% of restitution, contact time and max deflection.


clear data;
data.D = 50;
data.Quant = 100;
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
data.Bo = 0.017;
data.Oh = 0.006;
%Ang = 180;
%U = 18; %linspace(28, 50, 5)';

files = dir("**/simulation_postprocessing.mat");
currfol = pwd;
Westar = []; Bo = []; Oh = []; max_deflection = []; contact_time = []; 
coef_restitution = []; N = [];
plotting_data = table(Westar, Bo, Oh, max_deflection, contact_time, coef_restitution, N);
standalone = false;
saving = true;

for ii = 1:length(files)
    
    try 
        if length(files) == 1; folder_name = files.folder; else; folder_name = files(ii).folder; end
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
        load('muair.mat')
        load('g.mat','g') %gravitational constant
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
if ~isempty(plotting_data)
    if max(plotting_data.Bo) ~= min(plotting_data.Bo) || max(plotting_data.Oh) ~= min(plotting_data.Oh)
        warning("Bohn and Oh numbers not uniform across data");
    end

    close all
    cd(currfol);
    if standalone
    
        % Coef of Restitution
        figure(1);
        coef_res = scatter(plotting_data.Westar,plotting_data.coef_restitution,'MarkerEdgeColor',[ 0.4660    0.6740    0.1880],'LineWidth',4);
        grid on;
        %Center = plot(tvec(1:length(z)),z,'k','LineWidth',4);
        set(gca,'FontSize',16); %,'xlim',[0 16],'ylim',[-2 8])
        xlabel('   $We$   ','interpreter','LaTeX','FontSize',20); xlim([0, 2]);
        ylabel('   $\alpha \ \ \ $    ','interpreter','LaTeX','FontSize',20,'Rotation',0); ylim([0, 0.6]);
        title(sprintf("Coefficient of restitution for \n Bo = %.2g, Oh = %.2e", plotting_data.Bo(1), plotting_data.Oh(1)),'interpreter','LaTeX','FontSize',20);


        % Maximum deflection
        figure(2);
        max_def = scatter(plotting_data.Westar,plotting_data.max_deflection,'MarkerEdgeColor',[ 0.4660    0.6740    0.1880],'LineWidth',4);
        grid on;
        %Center = plot(tvec(1:length(z)),z,'k','LineWidth',4);
        set(gca,'FontSize',16); %,'xlim',[0 16],'ylim',[-2 8])
        xlabel('   $We$   ','interpreter','LaTeX','FontSize',20); xlim([0, 2]);
        ylabel('   $z / R_o \ \ \ $    ','interpreter','LaTeX','FontSize',20,'Rotation',0); ylim([0, 2.5]);
        title(sprintf("Maximum deflection for \n Bo = %.2g, Oh = %.2e", plotting_data.Bo(1), plotting_data.Oh(1)),'interpreter','LaTeX','FontSize',20);


        % ontact time
        figure(3);
        contact_time = scatter(plotting_data.Westar,plotting_data.contact_time,'MarkerEdgeColor',[ 0.4660    0.6740    0.1880],'LineWidth',4);
        grid on;
        %Center = plot(tvec(1:length(z)),z,'k','LineWidth',4);
        set(gca,'FontSize',16); %,'xlim',[0 16],'ylim',[-2 8])
        xlabel('   $We$   ','interpreter','LaTeX','FontSize',20); xlim([0, 2]);
        ylabel('   $t^* \ \ \ $    ','interpreter','LaTeX','FontSize',20,'Rotation',0); ylim([0, 6]);
        title(sprintf("Contact time for \n Bo = %.2g, Oh = %.2e", plotting_data.Bo(1), plotting_data.Oh(1)),'interpreter','LaTeX','FontSize',20);
    else
        
        if isfield(data, "Bo") && data.Bo > 0.03
            figname = "oil_3panel_FINALF.fig";
        else
            figname = "water_QPExp_DNS_3panel_FFF.fig";        
        end
        
        b = openfig(figname);
        b = b.Children;
        

        scatter(b(1), plotting_data.Westar,plotting_data.max_deflection,'MarkerEdgeColor',  [ 0.4660    0.6740    0.1880],'LineWidth',4);
        scatter(b(2), plotting_data.Westar,plotting_data.contact_time,'MarkerEdgeColor',    [ 0.4660    0.6740    0.1880],'LineWidth',4);
        scatter(b(3), plotting_data.Westar,plotting_data.coef_restitution,'MarkerEdgeColor',[ 0.4660    0.6740    0.1880],'LineWidth',4);
        
        if saving == true
            id = datetime('now'); id.Format = 'yyyyMMddmmss';

            f_maxdef = figure(2);
            copyobj(b(1), f_maxdef); a = gca;
            a.Position = [0.2, 0.1, 0.6, 0.9];
            saveas(f_maxdef, sprintf("../0_data/manual/maximum_deflection%s", id), 'eps');
            savefig(f_maxdef, sprintf("../0_data/manual/maximum_deflection%s.fig", id));

            f_contact = figure(3);
            copyobj(b(2), f_contact);a = gca;
            a.Position = [0.2, 0.1, 0.6, 0.9];
            saveas(f_contact, sprintf("../0_data/manual/contact_time%s", id), 'eps');
            savefig(f_contact, sprintf("../0_data/manual/contact_time%s.fig", id));

            f_CR = figure(4);
            copyobj(b(3), f_CR);a = gca;
            a.Position = [0.2, 0.1, 0.6, 0.9];
            saveas(f_CR, sprintf("../0_data/manual/coef_res%s", id), 'eps');
            savefig(f_CR, sprintf("../0_data/manual/coef_res%s.fig", id));
        end
    end
else
    warning("Couldn't find any simulations with the specified parameters");
end
function bool = is_valid(simul, data)
    bool = true; tolerance = 0.05;
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
