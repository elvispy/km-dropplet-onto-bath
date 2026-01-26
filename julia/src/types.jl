const NEAR_SINGULAR = Ref(false)
const NEAR_SINGULAR_TOL = 1e-12

struct CaseData
    Ro::Float64
    rhoS::Float64
    sigmaS::Float64
    rho::Float64
    sigma::Float64
    nu::Float64
    muair::Float64
    g::Float64
    D::Float64
    quant::Int
    nr::Int
    dr::Float64
    Delta::Matrix{Float64}
    IntMat::Matrix{Float64}
    DTN::Matrix{Float64}
    Ma::Float64
    Ra::Float64
end

mutable struct Condition
    deformation_amplitudes::Vector{Float64}
    deformation_velocities::Vector{Float64}
    pressure_amplitudes::Vector{Float64}
    dt::Float64
    nb_harmonics::Int
    current_time::Float64
    center_of_mass::Float64
    center_of_mass_velocity::Float64
    nb_contact_points::Float64
end

scalar(x) = x isa AbstractArray ? x[1] : x
eye(n::Int) = Matrix{Float64}(I, n, n)
