% This sript will try to sweep in the folder Structure all the empty
% folders

% List of files needed to perform the operation
aux_files = [ ...
    "zs_from_spherical.m", "r_from_spherical.m", ...
    "maximum_contact_radius.m", "theta_from_cylindrical.m", ...
    "project_amplitudes.m", "solve_ODE_unkown.m" , ...
    "calculate_tan.m", "solveDD0.m", "solvenDDCusp.m", ...
    "zeta_generator.m", "collectdnPl.m", "legendre_dx.m", ...
    "legendre_ddx.m", "collectPl.m", "my_legendre.m"];


% Initial folder to go back to
root = pwd;

% A folder which MUST have all the dependencies needed
safe_folder = fullfile(root, "D50Quant100\rho1000sigma7220nu98muair0\RhoS1000SigmaS7220\R0350mm\ImpDefCornerAng180U38");

% To force and repeat sweeps (.mat)
force_sweep = false;

files = dir("**/*ExactSH.m");
for ii = 1:length(files)
    cd(files(ii).folder);
    
    if force_sweep == true || isempty(dir("etaOri*.mat")) == true
        for file = aux_files
            if ~exist(file, "file")
                copyfile(fullfile(safe_folder, file), pwd)
            end
        end
        VertPolarExactSH;
    end
end
% file_names = {files.name};