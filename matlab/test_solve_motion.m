
function test_solve_motion()
    % TEST_SOLVE_MOTION Regression test for solve_motion.m
    %
    % This script runs both the original and the new solve_motion functions
    % on a small test case and asserts that the results (z, numl) are identical.
    
    % --- Configuration ---
    U0 = 10;
    N = 15;
    tolP = 1e-4; 
    
    % 1. Define paths
    % matlab/1_code/D25Quant100/rho870sigma1870nu0muair0.02/rhoS870sigmaS1870/R0250mm/ImpDefCornerAng180U10
    base_dir = fullfile(pwd, '1_code', 'D25Quant100', ...
        'rho870sigma1870nu0muair0.02', 'rhoS870sigmaS1870', 'R0250mm', ...
        'ImpDefCornerAng180U10');
        
    if ~exist(base_dir, 'dir')
        error('Base directory not found: %s', base_dir);
    end
    
    tolP_test = 1.00e-04;
    test_run_name = sprintf('N=%dtol=%0.2e', N, tolP_test);
    test_dir = fullfile(base_dir, test_run_name);
    
    % 2. Setup Test Directory
    if ~exist(test_dir, 'dir')
        mkdir(test_dir);
    end
    
    % Copy Ma.mat and Ra.mat from parent's parent if needed? 
    % load_simulation_setup looks in `solid_dir` for Ma/Ra.
    % solid_dir is parent of R0, i.e. RhoS...
    % It seems load_simulation_setup is fine as long as parents exist.
    % It only checks `exist(run_dir, 'dir')`.
    
    fprintf('Running tests in: %s\n', test_dir);
    
    % 3. Run Original
    fprintf('--- Running solve_motion_original ---\n');
    t_end_test = 0.5; % Run for a very short time (approx 10-20 steps usually if dt ~ 1e-3)
                      % Standard dt is small in dim units. T_unit ~ sqrt(rho R^3 / sigma)
                      % For water drop R=0.35mm:
                      % rho=1, R=0.035, sigma=72.
                      % T = sqrt(1 * 0.035^3 / 72) ~ sqrt(4e-5 / 72) ~ sqrt(0.5e-6) ~ 7e-4 s.
                      % tend=20 -> 20 * 0.7 ms = 14 ms.
                      % t_end_test = 0.2 -> 0.14 ms.
                      
    % We rename/backup existing files in test_dir if any
    clean_dir(test_dir);
    
    addpath(fullfile(pwd, '1_code', 'simulation_code'));
    
    % Run Original
    % We need to output to a specific place or rename files after.
    % solve_motion writes to `wd`.
    % We'll use a subdir for original and new?
    % load_simulation_setup REQURIES `wd` to be the `run_dir`. 
    % We cannot run in a subdir. 
    % We must run in `test_dir`, move files, then run again.
    
    try
        solve_motion_original(U0, [], N, tolP_test, test_dir, false, t_end_test);
    catch ME
        fprintf('Original crashed: %s\n', ME.message);
        rethrow(ME);
    end
    
    move_results(test_dir, 'original');
    
    % 4. Run New
    fprintf('--- Running solve_motion (new) ---\n');
    
    try
        solve_motion(U0, [], N, tolP_test, test_dir, true, t_end_test);
    catch ME
        fprintf('New version crashed: %s\n', ME.message);
        rethrow(ME);
    end
    
    move_results(test_dir, 'new');
    
    % 5. Compare
    compare_results(test_dir);
    
end

function clean_dir(d)
    files = dir(fullfile(d, '*.mat'));
    for i = 1:numel(files)
        delete(fullfile(d, files(i).name));
    end
end

function move_results(d, suffix)
    files = {'z.mat', 'numl.mat', 'oscillation_amplitudes.mat', 'pressure_amplitudes.mat'};
    for i = 1:numel(files)
        src = fullfile(d, files{i});
        if exist(src, 'file')
            dst = fullfile(d, [files{i}, '.', suffix]);
            movefile(src, dst);
        else
            warning('File %s not produced by %s run.', files{i}, suffix);
        end
    end
end

function compare_results(d)
    files = {'z.mat', 'numl.mat', 'oscillation_amplitudes.mat', 'pressure_amplitudes.mat'};
    failed = false;
    
    for i = 1:numel(files)
        f = files{i};
        f_orig = fullfile(d, [f, '.original']);
        f_new  = fullfile(d, [f, '.new']);
        
        if ~exist(f_orig, 'file') || ~exist(f_new, 'file')
             fprintf('SKIP: %s (missing files)\n', f);
             continue;
        end
        
        d_orig = load(f_orig);
        d_new = load(f_new);
        
        % extract the variable (same name as filename usually, minus .mat)
        [~, varname] = fileparts(f);
        if isfield(d_orig, varname)
            data_orig = d_orig.(varname);
            data_new = d_new.(varname);
            
            % Compare
            diff = norm(data_orig(:) - data_new(:));
            rel_diff = diff / (norm(data_orig(:)) + 1e-15);
            
            if rel_diff > 1e-10
                fprintf('FAIL: %s mismatch. Rel diff: %e\n', f, rel_diff);
                failed = true;
            else
                fprintf('PASS: %s matches. Rel diff: %e\n', f, rel_diff);
            end
        else
            warning('Variable %s not found in .mat file', varname);
        end
    end
    
    if failed
        error('Regression test FAILED.');
    else
        fprintf('Regression test PASSED.\n');
    end
end
