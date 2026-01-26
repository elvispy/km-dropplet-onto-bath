using DelimitedFiles
using JLD2

function read_csv(path::AbstractString)
    return readdlm(path, ',', Float64)
end

function read_scalar(path::AbstractString)
    return read_csv(path)[1]
end

function build_case_from_csv(csv_dir::AbstractString, out_path::AbstractString)
    data = Dict(
        "Ro" => read_scalar(joinpath(csv_dir, "Ro.csv")),
        "rhoS" => read_scalar(joinpath(csv_dir, "rhoS.csv")),
        "sigmaS" => read_scalar(joinpath(csv_dir, "sigmaS.csv")),
        "rho" => read_scalar(joinpath(csv_dir, "rho.csv")),
        "sigma" => read_scalar(joinpath(csv_dir, "sigma.csv")),
        "nu" => read_scalar(joinpath(csv_dir, "nu.csv")),
        "muair" => read_scalar(joinpath(csv_dir, "muair.csv")),
        "g" => read_scalar(joinpath(csv_dir, "g.csv")),
        "D" => read_scalar(joinpath(csv_dir, "D.csv")),
        "quant" => read_scalar(joinpath(csv_dir, "quant.csv")),
        "nr" => read_scalar(joinpath(csv_dir, "nr.csv")),
        "dr" => read_scalar(joinpath(csv_dir, "dr.csv")),
        "Delta" => read_csv(joinpath(csv_dir, "Delta.csv")),
        "IntMat" => read_csv(joinpath(csv_dir, "IntMat.csv")),
        "DTN" => read_csv(joinpath(csv_dir, "DTN.csv")),
        "Ma" => read_scalar(joinpath(csv_dir, "Ma.csv")),
        "Ra" => read_scalar(joinpath(csv_dir, "Ra.csv")),
    )

    JLD2.jldopen(out_path, "w") do file
        for (k, v) in data
            file[k] = v
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 2
        error("Usage: build_case_from_csv.jl <csv_dir> <out_path>")
    end
    build_case_from_csv(ARGS[1], ARGS[2])
end
