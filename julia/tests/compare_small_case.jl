using DelimitedFiles
using JLD2
using Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))

include(joinpath(ROOT, "src", "SolveMotion.jl"))

function normalize_array(val)
    if val isa AbstractVector
        return vec(val)
    end
    if val isa AbstractMatrix && (size(val, 1) == 1 || size(val, 2) == 1)
        return vec(val)
    end
    return val
end

function assert_close(name, julia_val, matlab_val; rtol=1e-8, atol=1e-10)
    jv = normalize_array(julia_val)
    mv = normalize_array(matlab_val)
    size(jv) == size(mv) || error("Size mismatch for $(name)")
    if !all(isapprox.(jv, mv; rtol=rtol, atol=atol))
        max_err = maximum(abs.(jv .- mv))
        error("Mismatch for $(name), max error $(max_err)")
    end
end

function read_csv(path::AbstractString)
    return readdlm(path, ',', Float64)
end

function escape_matlab_path(path::AbstractString)
    return replace(path, "'" => "''")
end

function run_test()
    U0 = parse(Float64, get(ENV, "TEST_U0", "10.0"))
    N = 10
    tolP = 1e-2
    max_steps = 5

    out_dir = mktempdir()
    matlab_out = joinpath(out_dir, "matlab")
    mkpath(matlab_out)

    matlab_sim = joinpath(ROOT, "..", "matlab", "1_code", "simulation_code")
    case_base = joinpath(ROOT, "..", "matlab", "1_code", "D5Quant20", "rho1000sigma7220nu98muair0",
        "RhoS1000SigmaS7220", "R0350mm")
    imp_dir = joinpath(case_base, @sprintf("ImpDefCornerAng180U%.4g", U0))
    run_dir = joinpath(imp_dir, @sprintf("N=%dtol=%0.2e", N, tolP))

    matlab_script = joinpath(out_dir, "matlab_run.m")
    open(matlab_script, "w") do io
        write(io, """
addpath('$(escape_matlab_path(matlab_sim))');
setenv('SOLVE_MOTION_MAX_STEPS','$(max_steps)');
setenv('SOLVE_MOTION_OUTPUT_PREFIX','matlab_');
if exist('$(escape_matlab_path(imp_dir))','dir') ~= 7, mkdir('$(escape_matlab_path(imp_dir))'); end
if exist('$(escape_matlab_path(run_dir))','dir') ~= 7, mkdir('$(escape_matlab_path(run_dir))'); end
delete(fullfile('$(escape_matlab_path(run_dir))', 'matlab_*'));
solve_motion($(U0), [], $(N), $(tolP), '$(escape_matlab_path(run_dir))', false);
out_dir = '$(escape_matlab_path(matlab_out))';
vars = {'z','etaOri','etas','vz','tvec','nlmax','numl','oscillation_amplitudes','pressure_amplitudes','Rv'};
for ii = 1:length(vars)
    name = vars{ii};
    f1 = fullfile('$(escape_matlab_path(run_dir))', ['matlab_' name '.mat']);
    f2 = fullfile('$(escape_matlab_path(run_dir))', ['matlab_.mat' name]);
    if exist(f1,'file')
        data = load(f1, name);
        delete(f1);
    else
        data = load('-mat', f2, name);
        delete(f2);
    end
    writematrix(data.(name), fullfile(out_dir, [name '.csv']));
end
delete(fullfile('$(escape_matlab_path(run_dir))', 'matlab_ProblemConditions.mat'));
""")
    end
    run(`matlab -nodisplay -nojvm -batch "run('$(escape_matlab_path(matlab_script))')"` )

    ENV["SOLVE_MOTION_OUTPUT_DIR"] = out_dir
    ENV["SOLVE_MOTION_OUTPUT_PREFIX"] = "julia_"
    ENV["SOLVE_MOTION_MAX_STEPS"] = string(max_steps)
    ENV["SOLVE_MOTION_RANDOM_SEED"] = "1234"

    SolveMotion.solve_motion(U0, nothing, N, tolP, nothing, false)

    vars = ["z","etaOri","etas","vz","tvec","nlmax","numl","oscillation_amplitudes","pressure_amplitudes","Rv"]
    sample_csv = joinpath(matlab_out, "z.csv")
    isfile(sample_csv) || error("MATLAB outputs missing: $(sample_csv)")

    julia_results = JLD2.load(joinpath(out_dir, "julia_results.jld2"))
    for name in vars
        matlab_val = read_csv(joinpath(matlab_out, name * ".csv"))
        julia_val = julia_results[name]
        assert_close(name, julia_val, matlab_val)
    end

    println("MATLAB and Julia results match.")
end

run_test()
