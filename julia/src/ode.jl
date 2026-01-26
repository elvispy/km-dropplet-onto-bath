function solve_ODE_unkown(deformation, pressures, dt, previous_conditions, PROBLEM_CONSTANTS)
    deformation_vec = deformation isa AbstractArray ? vec(deformation) : [deformation]
    pressures_vec = pressures isa AbstractArray ? vec(pressures) : [pressures]
    nb_harmonics = previous_conditions[end].nb_harmonics

    calc_vel = true
    if any(isnan.(deformation_vec)) || length(deformation_vec) < nb_harmonics
        deformation_vec = zeros(nb_harmonics)
        calc_vel = true
    end
    if any(isnan.(pressures_vec)) || length(pressures_vec) < nb_harmonics
        pressures_vec = zeros(nb_harmonics)
        calc_vel = false
    end

    n = length(previous_conditions)
    if n > 2 || n < 1
        error("Unexpected condition history size")
    end

    previous_velocities = zeros(Float64, nb_harmonics, n)
    previous_deformation = zeros(Float64, nb_harmonics, n)
    for ii in 1:n
        previous_velocities[:, ii] .= previous_conditions[ii].deformation_velocities
        previous_deformation[:, ii] .= previous_conditions[ii].deformation_amplitudes
    end

    if n == 1
        coefs = [-1.0, 1.0]
    else
        rk = dt / previous_conditions[end].dt
        ak = (1 + 2 * rk) / (1 + rk)
        bk = -(1 + rk)
        ck = rk^2 / (1 + rk)
        coefs = [ck, bk, ak]
    end

    harmonics = collect(1:nb_harmonics)
    amplitudes_coefficients = harmonics .* (harmonics .+ 2) .* (harmonics .- 1)
    pressure_coeffs = harmonics
    Oh = PROBLEM_CONSTANTS["Oh"]
    vel_coeffs = 2 .* (2 .* harmonics .+ 1) .* (harmonics .- 1)
    coef_end = coefs[end]
    coef_prev = reshape(coefs[1:(end - 1)], 1, :)
    result = (coef_end / dt .* (coef_end .+ dt .* Oh .* vel_coeffs) .+ amplitudes_coefficients .* dt) .* deformation_vec .+
        dt .* (pressure_coeffs .* pressures_vec) .+
        sum(coef_prev .* ((coef_end .+ dt .* Oh .* vel_coeffs) ./ dt .* previous_deformation .+ previous_velocities), dims=2)

    if calc_vel
        unkn = -vec(result) ./ (coef_end / dt .* (coef_end .+ dt .* Oh .* vel_coeffs) .+ amplitudes_coefficients .* dt)
        vel = (coef_end .* unkn .+ sum(coef_prev .* previous_deformation, dims=2)[:, 1]) ./ dt
        vel[1] = 0.0
        unkn[1] = 0.0
    else
        unkn = -vec(result) ./ (dt .* harmonics)
        vel = fill(NaN, 1)
    end
    return unkn, vel
end

