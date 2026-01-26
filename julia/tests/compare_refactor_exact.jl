using JLD2

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const RUN_SCRIPT = joinpath(ROOT, "julia", "tests", "run_solver_case.jl")

function arg_value(flag::String, default::String)
    idx = findfirst(==(flag), ARGS)
    return idx === nothing ? default : ARGS[idx + 1]
end

function ensure_baseline_tree(root::AbstractString)
    if haskey(ENV, "BASELINE_TREE")
        return ENV["BASELINE_TREE"]
    end
    commit = get(ENV, "BASELINE_COMMIT", "6441f10")
    worktree_dir = arg_value("--baseline-tree", mktempdir())
    run(`git -C $root worktree add $worktree_dir $commit`)
    return worktree_dir
end

function run_case(tree::AbstractString, out_dir::AbstractString, prefix::AbstractString, case_path::AbstractString)
    cmd = `julia --project=$(joinpath(ROOT, "julia")) $RUN_SCRIPT -- --tree $tree --out $out_dir --prefix $prefix --case-path $case_path`
    run(cmd)
end

function assert_exact(name, a, b)
    size(a) == size(b) || error("Size mismatch for $(name)")
    for idx in eachindex(a)
        isequal(a[idx], b[idx]) || error("Mismatch for $(name) at $(idx)")
    end
end

baseline_tree = ensure_baseline_tree(ROOT)
current_tree = ROOT
case_path = joinpath(ROOT, "julia", "data", "case_d5q20.jld2")

baseline_out = mktempdir()
current_out = mktempdir()

run_case(baseline_tree, baseline_out, "baseline_", case_path)
run_case(current_tree, current_out, "current_", case_path)

baseline_results = JLD2.load(joinpath(baseline_out, "baseline_results.jld2"))
current_results = JLD2.load(joinpath(current_out, "current_results.jld2"))

vars = ["z", "etaOri", "etas", "vz", "tvec", "nlmax", "numl", "oscillation_amplitudes", "pressure_amplitudes", "Rv"]
for name in vars
    assert_exact(name, current_results[name], baseline_results[name])
end

println("Exact match for baseline vs refactor outputs.")
