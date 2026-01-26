const ROOT = normpath(joinpath(@__DIR__, "..", ".."))

function arg_value(flag::String, default::String)
    idx = findfirst(==(flag), ARGS)
    return idx === nothing ? default : ARGS[idx + 1]
end

tree = arg_value("--tree", ROOT)
out_dir = arg_value("--out", mktempdir())
prefix = arg_value("--prefix", "run_")
u0 = parse(Float64, arg_value("--u0", "10.0"))
n = parse(Int, arg_value("--n", "10"))
tolp = parse(Float64, arg_value("--tolp", "1e-2"))
max_steps = arg_value("--max-steps", "5")
seed = arg_value("--seed", "1234")
case_path = arg_value("--case-path", "")

ENV["SOLVE_MOTION_OUTPUT_DIR"] = out_dir
ENV["SOLVE_MOTION_OUTPUT_PREFIX"] = prefix
ENV["SOLVE_MOTION_MAX_STEPS"] = max_steps
ENV["SOLVE_MOTION_RANDOM_SEED"] = seed
if !isempty(case_path)
    ENV["SOLVE_MOTION_CASE_PATH"] = case_path
end

include(joinpath(tree, "julia", "src", "SolveMotion.jl"))

SolveMotion.solve_motion(u0, nothing, n, tolp, nothing, false)
