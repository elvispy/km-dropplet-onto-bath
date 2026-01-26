module SolveMotion

using Dates
using LinearAlgebra
using Printf
using Random
using JLD2

export solve_motion, load_case, build_case, generate_dtn_new345, build_domain

include("types.jl")
include("case.jl")
include("interp.jl")
include("legendre.jl")
include("geometry.jl")
include("projection.jl")
include("linear.jl")
include("contact.jl")
include("ode.jl")
include("utils.jl")
include("domain.jl")
include("dtn.jl")
include("case_builder.jl")
include("solver.jl")

end
