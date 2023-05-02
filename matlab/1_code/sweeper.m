% This sript will try to sweep simulations according to two
% rules:
% 1) Parameters set in this sweep
% 2) All empty folders will be swept. 

% STEP 1: Define which simulations are to be run. 

D = 50;
Quant = 100;
rho = 1; % must multiply by x1000
sigma = 72.20; % must multiply by x100
muair = 0;
RhoS = 1; % must multiply by x1000
SigmaS = 72.20; % must multiply by x100
R = linspace(0.02, 0.05, 5)'; % must multiply by x10
Ang = 180;
U = linspace(28, 50, 5)';

[Didx, Quantidx, rhoidx, sigmaidx, muairidx, ...
    RhoSidx, SigmaSidx, Ridx, Angidx, Uidx] = ...
    ndgrid(1:length(D), 1:length(Quant), 1:length(rho), 1:length(sigma), ...
    1:length(muair), 1:length(RhoS), 1:length(SigmaS), ...
    1:length(R), 1:length(Ang), 1:length(U));

cartesian_product = [D(Didx, :), Quant(Quantidx, :), rho(rhoidx, :), ...
    sigma(sigmaidx, :), muair(muairidx, :), RhoS(RhoSidx, :), SigmaS(SigmaSidx, :), ...
    R(Ridx, :), Ang(Angidx, :), U(Uidx, :)];
% Turn simulations into table
if isempty(cartesian_product) == true; cartesian_product = double.empty(0, 10); end
simulations_cgs = array2table(cartesian_product, ...
    "VariableNames", ["D", "Quant", "rho", "sigma", "muair", "RhoS", ...
    "SigmaS", "R", "Ang", "U"]);
% Now you can manually add any simulations that you would like to run, such
% as:
simulations_cgs = [simulations_cgs; ...
    {50, 100, 1, 72.20, 0, 1, 72.20, 0.001, 180, 1}];


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
