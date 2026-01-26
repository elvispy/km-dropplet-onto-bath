using JLD2

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT, "julia", "src", "SolveMotion.jl"))

function assert_exact(name, a, b)
    size(a) == size(b) || error("Size mismatch for $(name)")
    for idx in eachindex(a)
        isequal(a[idx], b[idx]) || error("Mismatch for $(name) at $(idx)")
    end
end

case_path = get(ENV, "SOLVE_MOTION_CASE_PATH", joinpath(ROOT, "julia", "data", "case_d5q20.jld2"))
ENV["SOLVE_MOTION_CASE_PATH"] = case_path
ENV["SOLVE_MOTION_MAX_STEPS"] = get(ENV, "SOLVE_MOTION_MAX_STEPS", "5")
ENV["SOLVE_MOTION_RANDOM_SEED"] = get(ENV, "SOLVE_MOTION_RANDOM_SEED", "1234")

baseline_out = mktempdir()
current_out = mktempdir()

ENV["SOLVE_MOTION_OUTPUT_DIR"] = baseline_out
ENV["SOLVE_MOTION_OUTPUT_PREFIX"] = "baseline_"
SolveMotion.solve_motion_old(10.0, nothing, 10, 1e-2, nothing, false)

ENV["SOLVE_MOTION_OUTPUT_DIR"] = current_out
ENV["SOLVE_MOTION_OUTPUT_PREFIX"] = "current_"
SolveMotion.solve_motion(10.0, nothing, 10, 1e-2, nothing, false)

baseline_results = JLD2.load(joinpath(baseline_out, "baseline_results.jld2"))
current_results = JLD2.load(joinpath(current_out, "current_results.jld2"))

vars = ["z", "etaOri", "etas", "vz", "tvec", "nlmax", "numl", "oscillation_amplitudes", "pressure_amplitudes", "Rv"]
for name in vars
    assert_exact(name, current_results[name], baseline_results[name])
end

println("Exact match for solver_old vs solver outputs.")
