% This sript will try to sweep simulations according to two
% rules:
% 1) Parameters set in this sweep
% 2) All empty folders will be swept. 

% STEP 1: Define which simulations are to be run. 

% Idea: try to use a struct to add product of discrete sets and also
% manually addind simulation cases

% STEP 2: Recursively create folders and scripts needed
% Exmple of folder structure:
% D(50)Quant(100)
%   - DomainMaker.m
%   - ParRadDTNStops.m
%   -> rho(1000)sigma(7220)
%       - BathMaker.m
%       -> rhoS(1000)SigmaS(7220)
%       - DropFluidMaker.m    
%           -> R(0350)mm 
%               - RoMaker
%               -> ImpDefCornerAng(180)U(28)
%                   - VertPolarExactSH.m
%                   - And many others ...

% List of files needed to run the simulation and that do not need to change
% from simul to simul.
aux_files = [ ...
    "zs_from_spherical.m", "r_from_spherical.m", ...
    "maximum_contact_radius.m", "theta_from_cylindrical.m", ...
    "project_amplitudes.m", "solve_ODE_unkown.m" , ...
    "calculate_tan.m", "solveDD0.m", "solvenDDCusp.m", ...
    "zeta_generator.m", "collectdnPl.m", "legendre_dx.m", ...
    "legendre_ddx.m", "collectPl.m", "my_legendre.m"];


% Initial folder to go back to
root = pwd;

% A folder which MUST have all the dependencies needed (and all its
% parents, too. 
safe_folder = fullfile(root, "D50Quant100\rho1000sigma7220nu98muair0\RhoS1000SigmaS7220\R0350mm\ImpDefCornerAng180U38");

% To force and repeat sweeps (.mat)
force_sweep = false;


% STEP 3: Actually run the simulations. 
files = dir("**/*ExactSH.m");
for ii = 1:length(files)
    cd(files(ii).folder);
    
    % Check if etaOri exists (the center of the bath)
    if force_sweep == true || isempty(dir("etaOri*.mat")) == true
        for file = aux_files
            if ~exist(file, "file")
                copyfile(fullfile(safe_folder, file), pwd)
            end
        end
        VertPolarExactSH;
    end
end
