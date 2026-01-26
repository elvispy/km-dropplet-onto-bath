using Dates
using JLD2

# Simple sweep runner + CSV summary.
# Run with: julia --project=julia julia/scripts/sweep.jl

include(joinpath(@__DIR__, "..", "src", "SolveMotion.jl"))

function run_sweep_and_summarize(; outdir="julia/output", U0s=[5.0, 10.0, 15.0], N=10, tolP=1e-2)
    mkpath(outdir)
    ENV["SOLVE_MOTION_OUTPUT_DIR"] = outdir

    for U0 in U0s
        prefix = "U$(U0)_N$(N)_tol$(tolP)_" * Dates.format(now(), "yyyymmdd_HHMMSS") * "_"
        ENV["SOLVE_MOTION_OUTPUT_PREFIX"] = prefix
        SolveMotion.solve_motion(U0, nothing, N, tolP, nothing, false)
    end

    # Build a CSV summary of runs in the output folder.
    entries = String[]
    for f in readdir(outdir; join=true)
        if endswith(f, "problem_conditions.jld2")
            d = load(f)
            push!(entries, join([basename(f), d["U0"], d["N"], d["tend"], d["simul_time"]], ","))
        end
    end
    open(joinpath(outdir, "summary.csv"), "w") do io
        write(io, "file,U0,N,tend,simul_time\n")
        for line in entries
            write(io, line * "\n")
        end
    end
end

run_sweep_and_summarize()
