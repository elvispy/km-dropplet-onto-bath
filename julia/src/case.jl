function get_required(data::Dict{String, Any}, key::String)
    haskey(data, key) || error("Missing $(key) in case file")
    return data[key]
end

function default_case_path()
    candidates = [
        joinpath(@__DIR__, "data", "case_d5q20.jld2"),
        joinpath(@__DIR__, "..", "data", "case_d5q20.jld2"),
        joinpath(@__DIR__, "..", "..", "data", "case_d5q20.jld2"),
    ]
    for path in candidates
        if isfile(path)
            return normpath(path)
        end
    end
    return normpath(candidates[end])
end

function load_case(path::AbstractString)
    case_path = isdir(path) ? joinpath(path, "case.jld2") : path
    isfile(case_path) || error("Case file not found at $(case_path)")
    data = JLD2.load(case_path)
    Ro = float(scalar(get_required(data, "Ro")))
    rhoS = float(scalar(get_required(data, "rhoS")))
    sigmaS = float(scalar(get_required(data, "sigmaS")))
    rho = float(scalar(get_required(data, "rho")))
    sigma = float(scalar(get_required(data, "sigma")))
    nu = float(scalar(get_required(data, "nu")))
    muair = float(scalar(get_required(data, "muair")))
    g = float(scalar(get_required(data, "g")))
    D = float(scalar(get_required(data, "D")))
    quant = Int(round(scalar(get_required(data, "quant"))))
    nr = Int(round(scalar(get_required(data, "nr"))))
    dr = float(scalar(get_required(data, "dr")))
    Delta = Float64.(get_required(data, "Delta"))
    IntMat = Float64.(get_required(data, "IntMat"))
    DTN = Float64.(get_required(data, "DTN"))
    Ma = float(scalar(get_required(data, "Ma")))
    Ra = float(scalar(get_required(data, "Ra")))
    return CaseData(Ro, rhoS, sigmaS, rho, sigma, nu, muair, g, D, quant, nr, dr, Delta, IntMat, DTN, Ma, Ra)
end

function load_case()
    return load_case(default_case_path())
end
