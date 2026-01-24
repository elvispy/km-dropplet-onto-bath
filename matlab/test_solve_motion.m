% test_solve_motion.m
% Regression test for solve_motion.m
% Refactoring verification.

% Add simulation code to path
addpath(fullfile(pwd, '1_code/simulation_code'));

% Parameters (matches an existing folder structure in D5Quant20)
D = 5;
Quant = 20;
rho = 1.0;
sigma = 72.2;
nu = 0.0098;
muair = 0;
RhoS = 1.0;
SigmaS = 72.2;
Ro = 0.035;
Ang = 180;
U0 = 20;
N = 5;
tolP = 1e-5;

% Construct paths
base_dir = fullfile(pwd, '1_code', sprintf('D%dQuant%d', D, Quant));
fluid_dir = fullfile(base_dir, sprintf('rho%dsigma%dnu%dmuair%g', round(1000*rho), round(100*sigma), round(10000*nu), muair));
solid_dir = fullfile(fluid_dir, sprintf('RhoS%dSigmaS%d', round(1000*RhoS), round(100*SigmaS)));
r0_dir = fullfile(solid_dir, sprintf('R0%gmm', Ro*10000));
impact_dir = fullfile(r0_dir, sprintf('ImpDefCornerAng%gU%g', Ang, U0));

run_name = sprintf('N=%dtol=%0.2e', N, tolP);
run_dir = fullfile(impact_dir, run_name);

if ~exist(run_dir, 'dir')
    mkdir(run_dir);
end

% Run simulation
% We override tend to a small value to keep it short.
% T_unit depends on params, but roughly:
% T ~ sqrt(1 * 0.035^3 / 72.2) ~ sqrt(4e-5 / 70) ~ sqrt(6e-7) ~ 8e-4 s.
% 10 timesteps ~ 10 * dt. dt ~ 1/N^(1.5).
% Just run for a small generic time.
tend_override = 0.005; % dimensionless time?
% In solve_motion: tend = tend_override.
% Default tend is 20.
% We use debug_flag = true to see output.

fprintf('Running test simulation in %s\n', run_dir);
try
    solve_motion(U0, [], N, tolP, run_dir, true, tend_override);
    fprintf('Test PASSED\n');
catch ME
    fprintf('Test FAILED: %s\n', ME.message);
    fprintf('Stack trace:\n');
    for k = 1:length(ME.stack)
        fprintf('File: %s, Line: %d, Name: %s\n', ME.stack(k).file, ME.stack(k).line, ME.stack(k).name);
    end
end
