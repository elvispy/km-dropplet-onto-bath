% This sript will try to sweep simulations according to two
% rules:
% 1) Parameters set in this sweep
% 2) All empty folders will be swept. 

% STEP 1: Define which simulations are to be run. 


D = 5;
Quant = 20;

rho = 1; % must multiply by x1000
sigma = 72.20; % must multiply by x100
nu = 9.78E-3; % Multiply by x10000
muair = 0;
RhoS = 1; % must multiply by x1000
SigmaS = 72.20; % must multiply by x100
R = 0.035; % linspace(0.02, 0.05, 5)'; % must multiply by x10
Ang = 180;
U = 17.3; %[44.1, 47.1, 52.1]'; % linspace(58, 8, 11)';
modes = 25;

[Didx, Quantidx, rhoidx, sigmaidx, muairidx, nuidx, ...
    RhoSidx, SigmaSidx, Ridx, Angidx, Uidx, modesidx] = ...
    ndgrid(1:length(D), 1:length(Quant), 1:length(rho), 1:length(sigma), ...
    1:length(muair), 1:length(nu), 1:length(RhoS), 1:length(SigmaS), ...
    1:length(R), 1:length(Ang), 1:length(U), 1:length(modes));

cartesian_product = [D(Didx, :), Quant(Quantidx, :), rho(rhoidx, :), ...
    sigma(sigmaidx, :), muair(muairidx, :), nu(nuidx, :), RhoS(RhoSidx, :), ...
    SigmaS(SigmaSidx, :), R(Ridx, :), Ang(Angidx, :), U(Uidx, :), modes(modesidx, :)];

% Turn simulations into table
if isempty(cartesian_product) == true; cartesian_product = double.empty(0, 10); end
simulations_cgs = array2table(cartesian_product, ...
    "VariableNames", ["D", "Quant", "rho", "sigma", "muair", "nu", "RhoS", ...
    "SigmaS", "R", "Ang", "U", "modes"]);
% Now you can manually add any simulations that you would like to run, such
% as:
%  simulations_cgs = [simulations_cgs; ...
%      {50, 100, 1, 72.20, 0, 9.78E-3, 1, 72.20, 0.001, 180, 1}];
nrow = size(simulations_cgs, 1);
simulations_cgs.folder = repmat("", nrow, 1);

% STEP 2: Recursively create folders and scripts needed
% Exmple of folder structure:
% D(50)Quant(100)
%   - DomainMaker.m
%   - ParRadDTNStops.m
%   -> rho(1000)sigma(7220)nu(98)muair(0)
%       - BathMaker.m
%       -> rhoS(1000)SigmaS(7220)
%       - DropFluidMaker.m    
%           -> R(0350)mm 
%               - RoMaker
%               -> ImpDefCornerAng(180)U(28)
%                   - VertPolarExactSH.m
%                   - And many others ...

for ii = 1:height(simulations_cgs)
    simulations_cgs.folder(ii)  = create_folder_stucture(simulations_cgs(ii, :));
end

% List of files needed to run the simulation and that do not need to change
% from simul to simul.
aux_files = [ ...
    "zs_from_spherical.m", "r_from_spherical.m", ...
    "maximum_contact_radius.m", "theta_from_cylindrical.m", ...
    "custom_project_amplitudes.m", "solve_ODE_unkown.m" , ...
    "calculate_tan.m", "solveDD0.m", "solvenDDCusp.m", ...
    "zeta_generator.m", "collectdnPl.m", "collectPl.m"];

% To force and repeat sweeps (.mat)
force_sweep = false;

% STEP 3: Actually run the simulations. 

% Initial folder to go back to
root = pwd;
% A folder which MUST have all the dependencies needed (and all its
% parents, too. 
safe_folder = fullfile(root, "D50Quant100", "rho1000sigma7220nu98muair0", ...
    "RhoS1000SigmaS7220", "R0350mm", "ImpDefCornerAng180U38");

for ii = 1:height(simulations_cgs)
    %Check if etaOri exists (the center of the bath)
    cd(simulations_cgs.folder(ii));

    if force_sweep == true || isempty(dir("oscillation*.mat")) == true
        for file = aux_files
            try
                copyfile(fullfile(safe_folder, file), pwd)    
            catch ME
                if simulations_cgs.U(ii) ~= 38
                   throw("Unexpected error ocurred"); 
                end
            end
        end
        try
            callVertPolarExactSH();
        catch ME
           fprintf("Couldn't run simulation with the following parameters: \n Velcity: %g \n Modes: %g \n", ...
                simulations_cgs.U(ii), simulations_cgs.modes(ii)); 
        end
    end
    

end

% files = dir("**/*ExactSH.m");
% for ii = 1:length(files)
%     cd(files(ii).folder);
%     
%     % Check if etaOri exists (the center of the bath)
%     if force_sweep == true || isempty(dir("oscillation*.mat")) == true
%         for file = aux_files
%             if ~exist(file, "file")
%                 copyfile(fullfile(safe_folder, file), pwd)
%             end
%         end
% 
%         
%             VertPolarExactSH;
%     end
% end

function final_folder = create_folder_stucture(entry)
    base =  pwd;
    safe_folder = fullfile(base, "D50Quant100");
    
    % Defining folder structure
    physical_space = sprintf("D%gQuant%g", entry.D, entry.Quant);
    fluid_parameters = sprintf("rho%gsigma%gnu%.2gmuair%g", entry.rho*1000, entry.sigma*100, entry.nu*10000, entry.muair);
    sphere_parameters = sprintf("rhoS%gsigmaS%g", entry.RhoS*1000, entry.SigmaS*100);
    radius_folder = sprintf("R%04.4gmm", entry.R*10000);
    velocity_folder = sprintf("ImpDefCornerAng%gU%.4g", entry.Ang, entry.U);

    if ~isfolder(physical_space)
        mkdir(physical_space);
    end

    cd(physical_space);
    if ~isfile("zdrop.mat")
        copyfile(fullfile(safe_folder, "ParRadDTNStops.m"), pwd);
        
        s = fileread(fullfile(safe_folder, "DomainMaker.m"));
        s = regexprep(s, "D = [^\s]+;", sprintf("D = %g;", entry.D));
        s = regexprep(s, "quant = [^\s]+;", sprintf("quant = %g;", entry.Quant));
       
        writeID = fopen("DomainMaker.m", 'w+');
        fprintf(writeID, "%s", s);
        fclose(writeID);
        DomainMaker;
        if ~isfile("DTN*.mat")
            error("Integration matrix for fluid contact not found. Maually run ParRadDTNStops.m");
        end
    end


    safe_folder = fullfile(safe_folder, "rho1000sigma7220nu98muair0");
    if ~isfolder(fluid_parameters)
        mkdir(fluid_parameters);
    end
    cd(fluid_parameters);
    if ~isfile("muair.mat")
        s = fileread(fullfile(safe_folder, "BathMaker.m"));
        s = regexprep(s, "rho = [^\s]+;", sprintf("rho = %g;", entry.rho));
        s = regexprep(s, "sigma = [^\s]+;", sprintf("sigma = %g;", entry.sigma));
        s = regexprep(s, "nu = [^\s]+;", sprintf("nu = %.2e;", entry.nu));
        s = regexprep(s, "muair = [^\s]+;", sprintf("muair = %g;", entry.muair));
        
        writeID = fopen("BathMaker.m", 'w+');
        fprintf(writeID, "%s", s);
        fclose(writeID);
        BathMaker;
    end


    safe_folder = fullfile(safe_folder, "RhoS1000SigmaS7220");
    if ~isfolder(sphere_parameters)
        mkdir(sphere_parameters);
       %  copyfile(fullfile(safe_folder, "DropFluidMaker.m"), pwd);
    end
    cd(sphere_parameters);
    if ~isfile("sigmaS.mat")
        s = fileread(fullfile(safe_folder, "DropFluidMaker.m"));
        s = regexprep(s, "rhoS = [^\s]+;", sprintf("rhoS = %g;", entry.RhoS));
        s = regexprep(s, "sigmaS = [^\s]+;", sprintf("sigmaS = %g;", entry.SigmaS));
        
        writeID = fopen("DropFluidMaker.m", 'w+');
        fprintf(writeID, "%s", s);
        fclose(writeID);
        DropFluidMaker;
    end

    safe_folder = fullfile(safe_folder, "R0350mm");
    if ~isfolder(radius_folder)
        mkdir(radius_folder);
    end
    cd(radius_folder);
    if ~isfile("Ro.mat")
        s = sprintf("Ro = %d; save('Ro.mat','Ro') % Ball radius in cm. This will be our unit length", entry.R);
        writeID = fopen("RoMaker.m", "w+");
        fprintf(writeID, "%s", s);
        fclose(writeID);
        RoMaker;
    end

    safe_folder = fullfile(safe_folder, "ImpDefCornerAng180U38");
    if ~isfolder(velocity_folder)
        mkdir(velocity_folder);    
    end
    cd(velocity_folder);
    
    % Modify This file accordingly
    s = fileread(fullfile(safe_folder, "VertPolarExactSH.m"));
    s = regexprep(s, "U0 = 38;", sprintf("U0 = %g;", entry.U));
    s = regexprep(s, "N = \d+;", sprintf("N = %g;", entry.modes));
    writeID = fopen("VertPolarExactSH.m", 'w+');
    fprintf(writeID, "%s", s);
    fclose(writeID);
    final_folder = pwd;

    cd(base);
end

function callVertPolarExactSH()
    VertPolarExactSH;
end

function folder = return_folder(entry)
    folder = sprintf("D%gQuant%grho%gsigma%gnu%.2gmuair%grhoS%gsigmaS%gR%04.4gmmImpDefCornerAng%gU%2.3g", ...
                entry.D, entry.Quant, entry.rho*1000, entry.sigma*100, ...
                entry.nu*10000, entry.muair, entry.RhoS*1000, entry.SigmaS*100, ...
                entry.R*10000, entry.Ang, entry.U);
end
