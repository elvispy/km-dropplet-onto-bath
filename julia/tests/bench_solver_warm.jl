using Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "SolveMotion.jl"))

function arg_value(flag::String, default::String)
    idx = findfirst(==(flag), ARGS)
    return idx === nothing ? default : ARGS[idx + 1]
end

u0 = parse(Float64, arg_value("--u0", "10.0"))
n = parse(Int, arg_value("--n", "10"))
tolp = parse(Float64, arg_value("--tolp", "1e-2"))
max_steps = arg_value("--max-steps", "5")
seed = arg_value("--seed", "1234")
case_path = arg_value("--case-path", joinpath(ROOT, "data", "case_d5q20.jld2"))

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")
ENV["SOLVE_MOTION_MAX_STEPS"] = max_steps
ENV["SOLVE_MOTION_RANDOM_SEED"] = seed
ENV["SOLVE_MOTION_CASE_PATH"] = case_path

function run_case(fn, label)
    out_dir = mktempdir()
    ENV["SOLVE_MOTION_OUTPUT_DIR"] = out_dir
    ENV["SOLVE_MOTION_OUTPUT_PREFIX"] = label * "_"
    fn(u0, nothing, n, tolp, nothing, false)
    return out_dir
end

run_case(SolveMotion.solve_motion_old, "old_warm")
run_case(SolveMotion.solve_motion, "new_warm")

t_old = @elapsed run_case(SolveMotion.solve_motion_old, "old")
t_new = @elapsed run_case(SolveMotion.solve_motion, "new")

println(@sprintf("Warm timings (seconds): old=%0.3f new=%0.3f", t_old, t_new))
