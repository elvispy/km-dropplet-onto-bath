using Logging

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "SolveMotion.jl"))

# Ensure headless plotting backends (if loaded elsewhere) do not try to open a GUI.
ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")

precompile(SolveMotion.solve_motion, (Float64, Any, Int, Float64, Any, Bool))

ENV["SOLVE_MOTION_MAX_STEPS"] = get(ENV, "SOLVE_MOTION_MAX_STEPS", "5")
ENV["SOLVE_MOTION_RANDOM_SEED"] = get(ENV, "SOLVE_MOTION_RANDOM_SEED", "1234")
ENV["SOLVE_MOTION_OUTPUT_PREFIX"] = "precompile_"
ENV["SOLVE_MOTION_OUTPUT_DIR"] = mktempdir()

SolveMotion.solve_motion(10.0, nothing, 10, 1e-2, nothing, false)

@info "Precompile smoke run complete."
