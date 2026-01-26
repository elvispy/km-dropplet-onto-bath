using Printf

include(joinpath(@__DIR__, "..", "src", "SolveMotion.jl"))
using .SolveMotion

function env_float(key::AbstractString, default)
    val = get(ENV, key, "")
    return isempty(val) ? default : parse(Float64, val)
end

function env_int(key::AbstractString, default)
    val = get(ENV, key, "")
    return isempty(val) ? default : parse(Int, val)
end

D = env_float("CASE_D", 5.0)
quant = env_int("CASE_QUANT", 20)
rho = env_float("CASE_RHO", 1.0)
sigma = env_float("CASE_SIGMA", 72.2)
nu = env_float("CASE_NU", 9.78e-3)
muair = env_float("CASE_MUAIR", 0.0)
g = env_float("CASE_G", 981.0)
rhoS = env_float("CASE_RHOS", 1.0)
sigmaS = env_float("CASE_SIGMAS", 72.2)
Ro = env_float("CASE_RO", 0.035)
refp = env_int("CASE_REFP", 10)
num_batches = env_int("CASE_NUM_BATCHES", 20)
threaded = get(ENV, "CASE_THREADED", "true") != "false"

out_path = get(ENV, "CASE_OUTPUT", joinpath(@__DIR__, "..", "data", "case_d5q20.jld2"))

println(@sprintf("Building case D=%g quant=%d nr=ceil(D*quant/2) refp=%d", D, quant, refp))
SolveMotion.build_case(
    D=D,
    quant=quant,
    rho=rho,
    sigma=sigma,
    nu=nu,
    muair=muair,
    g=g,
    rhoS=rhoS,
    sigmaS=sigmaS,
    Ro=Ro,
    refp=refp,
    num_batches=num_batches,
    threaded=threaded,
    output_path=out_path,
)

println("Wrote case file to $(out_path)")
