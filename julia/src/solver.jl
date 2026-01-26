function solve_motion(U0, _unused, N, tolP, wd, debug_flag)
    seed = parse(Int, get(ENV, "SOLVE_MOTION_RANDOM_SEED", "1234"))
    Random.seed!(seed)
    tstart = time()
    Ang = 180.0
    case_path = get(ENV, "SOLVE_MOTION_CASE_PATH", "")
    if isempty(case_path)
        if wd === nothing || (wd isa AbstractString && isempty(wd))
            case = load_case()
            case_path = default_case_path()
        else
            case = load_case(string(wd))
            case_path = string(wd)
        end
    else
        case = load_case(case_path)
    end

    root = normpath(joinpath(case_path, "..", ".."))
    out_dir = get(ENV, "SOLVE_MOTION_OUTPUT_DIR", joinpath(root, "output"))
    mkpath(out_dir)
    prefix = get(ENV, "SOLVE_MOTION_OUTPUT_PREFIX", "")

    Ro = case.Ro
    rhoS = case.rhoS
    sigmaS = case.sigmaS
    rho = case.rho
    sigma = case.sigma
    nu = case.nu
    muair = case.muair
    g = case.g
    D = case.D
    quant = case.quant
    nr = case.nr
    dr = case.dr
    Delta = case.Delta
    IntMat = case.IntMat
    DTN = case.DTN
    Ma = case.Ma
    Ra = case.Ra

    L_unit = Ro
    M_unit = rhoS * L_unit^3
    T_unit = sqrt(rhoS * Ro^3 / sigmaS)
    V_unit = L_unit / T_unit

    Dr = rhoS / rho
    Re = L_unit^2 / (nu * T_unit)
    Fr = L_unit / (g * T_unit^2)
    We = rho * L_unit^3 / (sigma * T_unit^2)
    WeS = rhoS * Ro^3 / (sigmaS * T_unit^2)
    Westar = rhoS * U0^2 * Ro / sigmaS
    Oh = nu * sqrt(rhoS / (sigmaS * Ro))
    Cang = (Ang / 180) * pi
    tend = 20.0

    t = 0.0
    etao = zeros(Float64, nr)
    phio = zeros(Float64, nr)

    nsteps = ceil(Int, 20 / (2 * pi) * N^(3 / 2))
    dtb = 1.0 / nsteps
    steps = ceil(Int, (tend - t) / dtb)
    max_steps_env = get(ENV, "SOLVE_MOTION_MAX_STEPS", "")
    if !isempty(max_steps_env)
        steps = parse(Int, max_steps_env)
        tend = t + steps * dtb
    end

    if !isempty(max_steps_env)
        tvec = collect(t:dtb:(t + steps * dtb))
    else
        tvec = collect(t:dtb:(tend + 1))
    end
    tvecOri = tvec

    dt = tvec[2] - tvec[1]
    indexes_to_save = zeros(Int, steps + 1)
    current_to_save = 2
    indexes_to_save[1] = 1

    etaOri = zeros(Float64, steps + 1)
    z = zeros(Float64, steps + 1)
    vz = zeros(Float64, steps + 1)
    numl = zeros(Float64, steps + 1)
    oscillation_amplitudes = zeros(Float64, N, steps + 1)
    pressure_amplitudes = zeros(Float64, N, steps + 1)
    Rv = -ones(Float64, steps + 1)
    oscillation_velocities = zeros(Float64, N, steps + 1)
    nlmax = zeros(Float64, steps + 1)
    etas = zeros(Float64, nr, steps + 1)
    etas[:, 1] .= etao
    phis = zeros(Float64, nr, steps + 1)
    phis[:, 1] .= phio
    psMatPer = [zeros(Float64, quant + 1) for _ in 1:nsteps]
    zs = zeros(Float64, nr)

    omegas_frequencies = sqrt.(collect(1:N) .* (collect(1:N) .+ 2) .* (collect(1:N) .- 1) ./ WeS)
    amplitudes_old = oscillation_amplitudes[:, 1]
    amplitudes_velocities_old = oscillation_velocities[:, 1]
    B_l_ps_old = zeros(Float64, N)

    z[1] = -1 * zs_from_spherical(pi, oscillation_amplitudes[:, 1])
    vz[1] = -abs(U0 / V_unit)

    current_conditions = Condition(copy(amplitudes_old), copy(amplitudes_velocities_old), copy(B_l_ps_old), dt, N, 0.0, z[1], vz[1], 0.0)
    previous_conditions = [deepcopy(current_conditions), deepcopy(current_conditions)]
    previous_conditions[1].current_time = previous_conditions[2].current_time - dt
    previous_conditions[1].center_of_mass_velocity = previous_conditions[2].center_of_mass_velocity + dt / Fr
    previous_conditions[1].center_of_mass = previous_conditions[2].center_of_mass - previous_conditions[2].center_of_mass_velocity * dt

    gfun = (tval, idx) -> current_conditions.deformation_amplitudes[idx] * cos(omegas_frequencies[idx] * tval) +
        current_conditions.deformation_velocities[idx] / (omegas_frequencies[idx] + 1e-30) * sin(omegas_frequencies[idx] * tval)

    for idx in 1:N
        previous_conditions[1].deformation_amplitudes[idx] = gfun(-dt, idx)
        previous_conditions[1].deformation_velocities[idx] = (gfun(0, idx) - gfun(-2 * dt / 1000, idx)) / (2 * dt / 1000)
    end

    tentative_index = 0
    going_back = 0
    errortan = zeros(Float64, 5, steps + 1)
    ps_accepted = Float64[]
    ps_old = ps_accepted
    old_conditions = previous_conditions[1]

    PROBLEM_CONSTANTS = Dict(
        "froude_nb" => Fr,
        "weber_nb" => We,
        "nb_harmonics" => N,
        "density_ratio" => Dr,
        "omegas_frequencies" => omegas_frequencies,
        "spatial_tol" => dr,
        "initial_dt" => dtb,
        "DEBUG_FLAG" => debug_flag,
        "linear_on_theta" => true,
        "Ra" => Ra,
        "interpolation_number" => 10,
        "Oh" => Oh,
        "Westar" => Westar,
    )

    println(@sprintf("Starting simulation with case %s", case_path))

    exit_flag = false
    try
        while (t < tend) && !exit_flag
            NEAR_SINGULAR[] = false
            tentative_index += 1
            needed = tentative_index + 1
            if needed > length(z)
                extra = needed - length(z)
                resize!(etaOri, needed)
                resize!(z, needed)
                resize!(vz, needed)
                resize!(numl, needed)
                resize!(nlmax, needed)
                resize!(Rv, needed)
                resize!(indexes_to_save, needed)
                etaOri[(end - extra + 1):end] .= 0
                z[(end - extra + 1):end] .= 0
                vz[(end - extra + 1):end] .= 0
                numl[(end - extra + 1):end] .= 0
                nlmax[(end - extra + 1):end] .= 0
                Rv[(end - extra + 1):end] .= -1
                indexes_to_save[(end - extra + 1):end] .= 0
                etas = hcat(etas, zeros(Float64, nr, extra))
                phis = hcat(phis, zeros(Float64, nr, extra))
                errortan = hcat(errortan, 4 .* ones(Float64, 5, extra))
                oscillation_amplitudes = hcat(oscillation_amplitudes, zeros(Float64, N, extra))
                pressure_amplitudes = hcat(pressure_amplitudes, zeros(Float64, N, extra))
                oscillation_velocities = hcat(oscillation_velocities, zeros(Float64, N, extra))
            end
            t = tvec[tentative_index + 1]
            dt = t - tvec[tentative_index]

            if PROBLEM_CONSTANTS["DEBUG_FLAG"]
                println(@sprintf("Outside %0.4g, %0.3e", t - dt, dt))
            end

            etaprob = zeros(Float64, nr, 5)
            phiprob = zeros(Float64, nr, 5)
            vzprob = zeros(Float64, 5)
            zprob = zeros(Float64, 5)
            errortan[:, tentative_index + 1] .= 4

            psTent = ps_accepted
            RmaxOld = r_from_spherical(maximum_contact_radius(oscillation_amplitudes[:, tentative_index]),
                oscillation_amplitudes[:, tentative_index])
            nlmax[tentative_index] = max(floor(RmaxOld / dr) + 1, numl[tentative_index])
            nlmax_int = Int(round(nlmax[tentative_index]))

            thetaVec = theta_from_cylindrical(dr .* collect(0:(nlmax_int - 1)),
                oscillation_amplitudes[:, tentative_index])

            if norm(psTent, 1) == 0
                B_l_ps_tent = zeros(Float64, N)
            else
                nb_contact_points = Int(round(numl[tentative_index]))
                if PROBLEM_CONSTANTS["linear_on_theta"]
                    if nb_contact_points > 1
                        contactAngle = 1.5 * thetaVec[nb_contact_points] - 0.5 * thetaVec[nb_contact_points - 1]
                    else
                        contactAngle = (thetaVec[2] + thetaVec[1]) / 2
                    end
                    angles = range(contactAngle, stop=pi, length=(nb_contact_points + 1) * PROBLEM_CONSTANTS["interpolation_number"])
                    values = r_from_spherical(collect(angles), oscillation_amplitudes[:, tentative_index])
                    f_interp = r -> interp1_linear(dr .* collect(0:nb_contact_points),
                        vcat(psTent[1:nb_contact_points], 0.0), r, extrap=0.0)
                    values = f_interp.(values)
                    B_l_ps_tent = custom_project_amplitudes(collect(angles), values, N, NaN, NaN)
                else
                    error("project_amplitudes path is not implemented in this translation")
                end
            end

            amplitudes_tent, _ = solve_ODE_unkown(NaN, B_l_ps_tent, dt, previous_conditions, PROBLEM_CONSTANTS)

            RmaxTent = r_from_spherical(maximum_contact_radius(oscillation_amplitudes[:, tentative_index]),
                oscillation_amplitudes[:, tentative_index])
            nlmaxTent = floor(RmaxTent / dr) + 1
            nlmaxTent_int = Int(round(nlmaxTent))
            thetaVec = theta_from_cylindrical(dr .* collect(0:(nlmaxTent_int - 1)),
                oscillation_amplitudes[:, tentative_index])

            RvTent = zs_from_spherical(pi, amplitudes_tent)
            zs[1:nlmaxTent_int] .= zs_from_spherical(thetaVec, amplitudes_tent) .- RvTent
            zs[(nlmaxTent_int + 1):nr] .= Inf

            tanDrop = calculate_tan(dr .* collect(1:nlmaxTent_int) .- dr / 2, amplitudes_tent)
            angleDropMP = zeros(Float64, nlmaxTent_int)
            angleDropMP[1:nlmaxTent_int] .= atan.(tanDrop[1:nlmaxTent_int])

            psprob = zeros(Float64, nlmaxTent_int, 5)
            errorP = 1.0
            reduc = 0
            ll = 0

            while abs(errorP) >= tolP && reduc == 0
                ll += 1
                numl_curr = numl[tentative_index]

                if numl_curr < 0.5
                    etaprob[:, 3], phiprob[:, 3], zprob[3], vzprob[3], errortan[3, tentative_index + 1] =
                        solveDD0(dt, z[tentative_index], vz[tentative_index], etas[:, tentative_index],
                            phis[:, tentative_index], nr, Re, Delta, DTN, Fr, We, zs, RvTent)
                    if abs(errortan[3, tentative_index + 1]) < 0.5
                        numlTent = 0.0
                        etaTent = etaprob[:, 3]
                        phiTent = phiprob[:, 3]
                        psNew = zeros(Float64, nlmaxTent_int)
                        zTent = zprob[3]
                        vzTent = vzprob[3]
                    else
                        idx_prev = prev_index_for(numl, tentative_index, 1.0)
                        etaprob[:, 4], phiprob[:, 4], zprob[4], vzprob[4], ps_tmp, errortan[4, tentative_index + 1] =
                            solvenDDCusp(numl[idx_prev], 1, dt, z[tentative_index], vz[tentative_index],
                                etas[:, tentative_index], phis[:, tentative_index], nr, dr, Re, Delta, DTN, Fr,
                                We, Ma, zs, IntMat[1, :], angleDropMP, Cang, Dr, RvTent)
                        psprob[1, 4] = ps_tmp[1]
                        idx_prev = prev_index_for(numl, tentative_index, 2.0)
                        _, _, _, _, _, errortan[5, tentative_index + 1] =
                            solvenDDCusp(numl[idx_prev], 2, dt, z[tentative_index], vz[tentative_index],
                                etas[:, tentative_index], phis[:, tentative_index], nr, dr, Re, Delta, DTN, Fr,
                                We, Ma, zs, IntMat[2, :], angleDropMP, Cang, Dr, RvTent)
                        if abs(errortan[4, tentative_index + 1]) < abs(errortan[5, tentative_index + 1])
                            numlTent = 1.0
                            etaTent = etaprob[:, 4]
                            phiTent = phiprob[:, 4]
                            psNew = psprob[:, 4]
                            zTent = zprob[4]
                            vzTent = vzprob[4]
                        else
                            tvec = vcat(tvec[1:tentative_index],
                                (tvec[tentative_index] + tvec[tentative_index + 1]) / 2,
                                tvec[(tentative_index + 1):end])
                            tentative_index -= 1
                            reduc = 1
                        end
                    end
                elseif numl_curr > 0.5 && numl_curr < 1.5
                    etaprob[:, 2], phiprob[:, 2], zprob[2], vzprob[2], errortan[2, tentative_index + 1] =
                        solveDD0(dt, z[tentative_index], vz[tentative_index], etas[:, tentative_index],
                            phis[:, tentative_index], nr, Re, Delta, DTN, Fr, We, zs, RvTent)
                    if abs(errortan[2, tentative_index + 1]) < 0.5
                        numlTent = 0.0
                        etaTent = etaprob[:, 2]
                        phiTent = phiprob[:, 2]
                        psNew = zeros(Float64, nlmaxTent_int)
                        zTent = zprob[2]
                        vzTent = vzprob[2]
                    else
                        idx_prev = prev_index_for(numl, tentative_index, 1.0)
                        etaprob[:, 3], phiprob[:, 3], zprob[3], vzprob[3], ps_tmp, errortan[3, tentative_index + 1] =
                            solvenDDCusp(numl[idx_prev], 1, dt, z[tentative_index], vz[tentative_index],
                                etas[:, tentative_index], phis[:, tentative_index], nr, dr, Re, Delta, DTN, Fr,
                                We, Ma, zs, IntMat[1, :], angleDropMP, Cang, Dr, RvTent)
                        psprob[1, 3] = ps_tmp[1]
                        idx_prev = prev_index_for(numl, tentative_index, 2.0)
                        etaprob[:, 4], phiprob[:, 4], zprob[4], vzprob[4], psprob[1:2, 4], errortan[4, tentative_index + 1] =
                            solvenDDCusp(numl[idx_prev], 2, dt, z[tentative_index], vz[tentative_index],
                                etas[:, tentative_index], phis[:, tentative_index], nr, dr, Re, Delta, DTN, Fr,
                                We, Ma, zs, IntMat[2, :], angleDropMP, Cang, Dr, RvTent)
                        if abs(errortan[3, tentative_index + 1]) < abs(errortan[4, tentative_index + 1])
                            numlTent = 1.0
                            etaTent = etaprob[:, 3]
                            phiTent = phiprob[:, 3]
                            psNew = psprob[:, 3]
                            zTent = zprob[3]
                            vzTent = vzprob[3]
                        else
                            idx_prev = prev_index_for(numl, tentative_index, 3.0)
                            _, _, _, _, _, errortan[5, tentative_index + 1] =
                                solvenDDCusp(numl[idx_prev], 3, dt, z[tentative_index], vz[tentative_index],
                                    etas[:, tentative_index], phis[:, tentative_index], nr, dr, Re, Delta, DTN, Fr,
                                    We, Ma, zs, IntMat[3, :], angleDropMP, Cang, Dr, RvTent)
                            if abs(errortan[4, tentative_index + 1]) < abs(errortan[5, tentative_index + 1])
                                numlTent = 2.0
                                etaTent = etaprob[:, 4]
                                phiTent = phiprob[:, 4]
                                psNew = psprob[:, 4]
                                zTent = zprob[4]
                                vzTent = vzprob[4]
                            else
                                tvec = vcat(tvec[1:tentative_index],
                                    (tvec[tentative_index] + tvec[tentative_index + 1]) / 2,
                                    tvec[(tentative_index + 1):end])
                                tentative_index -= 1
                                reduc = 1
                            end
                        end
                    end
                elseif numl_curr > 1.5 && numl_curr < 2.5
                    _, _, _, _, errortan[1, tentative_index + 1] =
                        solveDD0(dt, z[tentative_index], vz[tentative_index], etas[:, tentative_index],
                            phis[:, tentative_index], nr, Re, Delta, DTN, Fr, We, zs, RvTent)
                    if abs(errortan[1, tentative_index + 1]) < 0.5
                        tvec = vcat(tvec[1:tentative_index],
                            (tvec[tentative_index] + tvec[tentative_index + 1]) / 2,
                            tvec[(tentative_index + 1):end])
                        tentative_index -= 1
                        reduc = 1
                    else
                        idx_prev = prev_index_for(numl, tentative_index, 2.0)
                        etaprob[:, 3], phiprob[:, 3], zprob[3], vzprob[3], psprob[1:2, 3], errortan[3, tentative_index + 1] =
                            solvenDDCusp(numl[idx_prev], 2, dt, z[tentative_index], vz[tentative_index],
                                etas[:, tentative_index], phis[:, tentative_index], nr, dr, Re, Delta, DTN, Fr,
                                We, Ma, zs, IntMat[2, :], angleDropMP, Cang, Dr, RvTent)
                        idx_prev = prev_index_for(numl, tentative_index, 1.0)
                        etaprob[:, 2], phiprob[:, 2], zprob[2], vzprob[2], ps_tmp, errortan[2, tentative_index + 1] =
                            solvenDDCusp(numl[idx_prev], 1, dt, z[tentative_index], vz[tentative_index],
                                etas[:, tentative_index], phis[:, tentative_index], nr, dr, Re, Delta, DTN, Fr,
                                We, Ma, zs, IntMat[1, :], angleDropMP, Cang, Dr, RvTent)
                        psprob[1, 2] = ps_tmp[1]
                        if abs(errortan[2, tentative_index + 1]) < abs(errortan[3, tentative_index + 1])
                            numlTent = 1.0
                            etaTent = etaprob[:, 2]
                            phiTent = phiprob[:, 2]
                            psNew = psprob[:, 2]
                            zTent = zprob[2]
                            vzTent = vzprob[2]
                        else
                            idx_prev = prev_index_for(numl, tentative_index, 3.0)
                            etaprob[:, 4], phiprob[:, 4], zprob[4], vzprob[4], psprob[1:3, 4], errortan[4, tentative_index + 1] =
                                solvenDDCusp(numl[idx_prev], 3, dt, z[tentative_index], vz[tentative_index],
                                    etas[:, tentative_index], phis[:, tentative_index], nr, dr, Re, Delta, DTN, Fr,
                                    We, Ma, zs, IntMat[3, :], angleDropMP, Cang, Dr, RvTent)
                            if abs(errortan[3, tentative_index + 1]) < abs(errortan[4, tentative_index + 1])
                                numlTent = 2.0
                                etaTent = etaprob[:, 3]
                                phiTent = phiprob[:, 3]
                                psNew = psprob[:, 3]
                                zTent = zprob[3]
                                vzTent = vzprob[3]
                            else
                                idx_prev = prev_index_for(numl, tentative_index, 4.0)
                                _, _, _, _, _, errortan[5, tentative_index + 1] =
                                    solvenDDCusp(numl[idx_prev], 4, dt, z[tentative_index], vz[tentative_index],
                                        etas[:, tentative_index], phis[:, tentative_index], nr, dr, Re, Delta, DTN, Fr,
                                        We, Ma, zs, IntMat[4, :], angleDropMP, Cang, Dr, RvTent)
                                if abs(errortan[4, tentative_index + 1]) < abs(errortan[5, tentative_index + 1])
                                    numlTent = 3.0
                                    etaTent = etaprob[:, 4]
                                    phiTent = phiprob[:, 4]
                                    psNew = psprob[:, 4]
                                    zTent = zprob[4]
                                    vzTent = vzprob[4]
                                else
                                    tvec = vcat(tvec[1:tentative_index],
                                        (tvec[tentative_index] + tvec[tentative_index + 1]) / 2,
                                        tvec[(tentative_index + 1):end])
                                    tentative_index -= 1
                                    reduc = 1
                                end
                            end
                        end
                    end
                elseif numl_curr > 2.5 && numl_curr < nlmaxTent - 1.5
                    idx_prev = prev_index_for(numl, tentative_index, numl_curr)
                    numl_int = Int(round(numl_curr))
                    etaprob[:, 3], phiprob[:, 3], zprob[3], vzprob[3], psprob[1:numl_int, 3], errortan[3, tentative_index + 1] =
                        solvenDDCusp(numl[idx_prev], numl_int, dt, z[tentative_index], vz[tentative_index],
                            etas[:, tentative_index], phis[:, tentative_index], nr, dr, Re, Delta, DTN, Fr,
                            We, Ma, zs, IntMat[numl_int, :], angleDropMP, Cang, Dr, RvTent)
                    idx_prev = prev_index_for(numl, tentative_index, numl_curr - 1)
                    etaprob[:, 2], phiprob[:, 2], zprob[2], vzprob[2], psprob[1:(numl_int - 1), 2], errortan[2, tentative_index + 1] =
                        solvenDDCusp(numl[idx_prev], numl_int - 1, dt, z[tentative_index], vz[tentative_index],
                            etas[:, tentative_index], phis[:, tentative_index], nr, dr, Re, Delta, DTN, Fr,
                            We, Ma, zs, IntMat[numl_int - 1, :], angleDropMP, Cang, Dr, RvTent)
                    if abs(errortan[2, tentative_index + 1]) < abs(errortan[3, tentative_index + 1])
                        idx_prev = prev_index_for(numl, tentative_index, numl_curr - 2)
                        _, _, _, _, _, errortan[1, tentative_index + 1] =
                            solvenDDCusp(numl[idx_prev], numl_int - 2, dt, z[tentative_index], vz[tentative_index],
                                etas[:, tentative_index], phis[:, tentative_index], nr, dr, Re, Delta, DTN, Fr,
                                We, Ma, zs, IntMat[numl_int - 2, :], angleDropMP, Cang, Dr, RvTent)
                        if abs(errortan[2, tentative_index + 1]) < abs(errortan[1, tentative_index + 1])
                            numlTent = numl_curr - 1
                            etaTent = etaprob[:, 2]
                            phiTent = phiprob[:, 2]
                            psNew = psprob[:, 2]
                            zTent = zprob[2]
                            vzTent = vzprob[2]
                        else
                            tvec = vcat(tvec[1:tentative_index],
                                (tvec[tentative_index] + tvec[tentative_index + 1]) / 2,
                                tvec[(tentative_index + 1):end])
                            tentative_index -= 1
                            reduc = 1
                        end
                    else
                        idx_prev = prev_index_for(numl, tentative_index, numl_curr + 1)
                        etaprob[:, 4], phiprob[:, 4], zprob[4], vzprob[4], psprob[1:(numl_int + 1), 4], errortan[4, tentative_index + 1] =
                            solvenDDCusp(numl[idx_prev], numl_int + 1, dt, z[tentative_index], vz[tentative_index],
                                etas[:, tentative_index], phis[:, tentative_index], nr, dr, Re, Delta, DTN, Fr,
                                We, Ma, zs, IntMat[numl_int + 1, :], angleDropMP, Cang, Dr, RvTent)
                        if abs(errortan[3, tentative_index + 1]) < abs(errortan[4, tentative_index + 1])
                            numlTent = numl_curr
                            etaTent = etaprob[:, 3]
                            phiTent = phiprob[:, 3]
                            psNew = psprob[:, 3]
                            zTent = zprob[3]
                            vzTent = vzprob[3]
                        else
                            idx_prev = prev_index_for(numl, tentative_index, numl_curr + 2)
                            _, _, _, _, _, errortan[5, tentative_index + 1] =
                                solvenDDCusp(numl[idx_prev], numl_int + 2, dt, z[tentative_index], vz[tentative_index],
                                    etas[:, tentative_index], phis[:, tentative_index], nr, dr, Re, Delta, DTN, Fr,
                                    We, Ma, zs, IntMat[numl_int + 2, :], angleDropMP, Cang, Dr, RvTent)
                            if abs(errortan[4, tentative_index + 1]) < abs(errortan[5, tentative_index + 1])
                                numlTent = numl_curr + 1
                                etaTent = etaprob[:, 4]
                                phiTent = phiprob[:, 4]
                                psNew = psprob[:, 4]
                                zTent = zprob[4]
                                vzTent = vzprob[4]
                            else
                                tvec = vcat(tvec[1:tentative_index],
                                    (tvec[tentative_index] + tvec[tentative_index + 1]) / 2,
                                    tvec[(tentative_index + 1):end])
                                tentative_index -= 1
                                reduc = 1
                            end
                        end
                    end
                elseif numl_curr > nlmax[tentative_index] - 1.5 && numl_curr < nlmaxTent - 0.5
                    idx_prev = prev_index_for(numl, tentative_index, numl_curr)
                    numl_int = Int(round(numl_curr))
                    etaprob[:, 3], phiprob[:, 3], zprob[3], vzprob[3], psprob[1:numl_int, 3], errortan[3, tentative_index + 1] =
                        solvenDDCusp(numl[idx_prev], numl_int, dt, z[tentative_index], vz[tentative_index],
                            etas[:, tentative_index], phis[:, tentative_index], nr, dr, Re, Delta, DTN, Fr,
                            We, Ma, zs, IntMat[numl_int, :], angleDropMP, Cang, Dr, RvTent)
                    idx_prev = prev_index_for(numl, tentative_index, numl_curr - 1)
                    etaprob[:, 2], phiprob[:, 2], zprob[2], vzprob[2], psprob[1:(numl_int - 1), 2], errortan[2, tentative_index + 1] =
                        solvenDDCusp(numl[idx_prev], numl_int - 1, dt, z[tentative_index], vz[tentative_index],
                            etas[:, tentative_index], phis[:, tentative_index], nr, dr, Re, Delta, DTN, Fr,
                            We, Ma, zs, IntMat[numl_int - 1, :], angleDropMP, Cang, Dr, RvTent)
                    if abs(errortan[2, tentative_index + 1]) < abs(errortan[3, tentative_index + 1])
                        idx_prev = prev_index_for(numl, tentative_index, numl_curr - 2)
                        _, _, _, _, _, errortan[1, tentative_index + 1] =
                            solvenDDCusp(numl[idx_prev], numl_int - 2, dt, z[tentative_index], vz[tentative_index],
                                etas[:, tentative_index], phis[:, tentative_index], nr, dr, Re, Delta, DTN, Fr,
                                We, Ma, zs, IntMat[numl_int - 2, :], angleDropMP, Cang, Dr, RvTent)
                        if abs(errortan[2, tentative_index + 1]) < abs(errortan[1, tentative_index + 1])
                            numlTent = numl_curr - 1
                            etaTent = etaprob[:, 2]
                            phiTent = phiprob[:, 2]
                            psNew = psprob[:, 2]
                            zTent = zprob[2]
                            vzTent = vzprob[2]
                        else
                            tvec = vcat(tvec[1:tentative_index],
                                (tvec[tentative_index] + tvec[tentative_index + 1]) / 2,
                                tvec[(tentative_index + 1):end])
                            tentative_index -= 1
                            reduc = 1
                        end
                    else
                        idx_prev = prev_index_for(numl, tentative_index, numl_curr + 1)
                        etaprob[:, 4], phiprob[:, 4], zprob[4], vzprob[4], psprob[1:(numl_int + 1), 4], errortan[4, tentative_index + 1] =
                            solvenDDCusp(numl[idx_prev], numl_int + 1, dt, z[tentative_index], vz[tentative_index],
                                etas[:, tentative_index], phis[:, tentative_index], nr, dr, Re, Delta, DTN, Fr,
                                We, Ma, zs, IntMat[numl_int + 1, :], angleDropMP, Cang, Dr, RvTent)
                        if abs(errortan[3, tentative_index + 1]) < abs(errortan[4, tentative_index + 1])
                            numlTent = numl_curr
                            etaTent = etaprob[:, 3]
                            phiTent = phiprob[:, 3]
                            psNew = psprob[:, 3]
                            zTent = zprob[3]
                            vzTent = vzprob[3]
                        else
                            numlTent = numl_curr + 1
                            etaTent = etaprob[:, 4]
                            phiTent = phiprob[:, 4]
                            psNew = psprob[:, 4]
                            zTent = zprob[4]
                            vzTent = vzprob[4]
                        end
                    end
                elseif numl_curr == nlmaxTent
                    idx_prev = prev_index_for(numl, tentative_index, numl_curr)
                    numl_int = Int(round(numl_curr))
                    etaprob[:, 3], phiprob[:, 3], zprob[3], vzprob[3], psprob[1:numl_int, 3], errortan[3, tentative_index + 1] =
                        solvenDDCusp(numl[idx_prev], numl_int, dt, z[tentative_index], vz[tentative_index],
                            etas[:, tentative_index], phis[:, tentative_index], nr, dr, Re, Delta, DTN, Fr,
                            We, Ma, zs, IntMat[numl_int, :], angleDropMP, Cang, Dr, RvTent)
                    idx_prev = prev_index_for(numl, tentative_index, numl_curr - 1)
                    etaprob[:, 2], phiprob[:, 2], zprob[2], vzprob[2], psprob[1:(numl_int - 1), 2], errortan[2, tentative_index + 1] =
                        solvenDDCusp(numl[idx_prev], numl_int - 1, dt, z[tentative_index], vz[tentative_index],
                            etas[:, tentative_index], phis[:, tentative_index], nr, dr, Re, Delta, DTN, Fr,
                            We, Ma, zs, IntMat[numl_int - 1, :], angleDropMP, Cang, Dr, RvTent)
                    if abs(errortan[2, tentative_index + 1]) < abs(errortan[3, tentative_index + 1])
                        idx_prev = prev_index_for(numl, tentative_index, numl_curr - 2)
                        _, _, _, _, _, errortan[1, tentative_index + 1] =
                            solvenDDCusp(numl[idx_prev], numl_int - 2, dt, z[tentative_index], vz[tentative_index],
                                etas[:, tentative_index], phis[:, tentative_index], nr, dr, Re, Delta, DTN, Fr,
                                We, Ma, zs, IntMat[numl_int - 2, :], angleDropMP, Cang, Dr, RvTent)
                        if abs(errortan[2, tentative_index + 1]) < abs(errortan[1, tentative_index + 1])
                            numlTent = numl_curr - 1
                            etaTent = etaprob[:, 2]
                            phiTent = phiprob[:, 2]
                            psNew = psprob[:, 2]
                            zTent = zprob[2]
                            vzTent = vzprob[2]
                        else
                            tvec = vcat(tvec[1:tentative_index],
                                (tvec[tentative_index] + tvec[tentative_index + 1]) / 2,
                                tvec[(tentative_index + 1):end])
                            tentative_index -= 1
                            reduc = 1
                        end
                    else
                        numlTent = numl_curr
                        etaTent = etaprob[:, 3]
                        phiTent = phiprob[:, 3]
                        psNew = psprob[:, 3]
                        zTent = zprob[3]
                        vzTent = vzprob[3]
                    end
                else
                    idx_prev = prev_index_for(numl, tentative_index, numl_curr - 1)
                    numl_int = Int(round(numl_curr))
                    etaprob[:, 2], phiprob[:, 2], zprob[2], vzprob[2], psprob[1:(numl_int - 1), 2], errortan[2, tentative_index + 1] =
                        solvenDDCusp(numl[idx_prev], numl_int - 1, dt, z[tentative_index], vz[tentative_index],
                            etas[:, tentative_index], phis[:, tentative_index], nr, dr, Re, Delta, DTN, Fr,
                            We, Ma, zs, IntMat[numl_int - 1, :], angleDropMP, Cang, Dr, RvTent)
                    if abs(errortan[2, tentative_index + 1]) < 4
                        idx_prev = prev_index_for(numl, tentative_index, numl_curr - 2)
                        _, _, _, _, _, errortan[1, tentative_index + 1] =
                            solvenDDCusp(numl[idx_prev], numl_int - 2, dt, z[tentative_index], vz[tentative_index],
                                etas[:, tentative_index], phis[:, tentative_index], nr, dr, Re, Delta, DTN, Fr,
                                We, Ma, zs, IntMat[numl_int - 2, :], angleDropMP, Cang, Dr, RvTent)
                        if abs(errortan[2, tentative_index + 1]) < abs(errortan[1, tentative_index + 1])
                            numlTent = numl_curr - 1
                            etaTent = etaprob[:, 2]
                            phiTent = phiprob[:, 2]
                            psNew = psprob[:, 2]
                            zTent = zprob[2]
                            vzTent = vzprob[2]
                        else
                            tvec = vcat(tvec[1:tentative_index],
                                (tvec[tentative_index] + tvec[tentative_index + 1]) / 2,
                                tvec[(tentative_index + 1):end])
                            tentative_index -= 1
                            reduc = 1
                        end
                    else
                        tvec = vcat(tvec[1:tentative_index],
                            (tvec[tentative_index] + tvec[tentative_index + 1]) / 2,
                            tvec[(tentative_index + 1):end])
                        tentative_index -= 1
                        reduc = 1
                    end
                end

                if ll == 100 && reduc == 0
                    tvec = vcat(tvec[1:tentative_index],
                        (tvec[tentative_index] + tvec[tentative_index + 1]) / 2,
                        tvec[(tentative_index + 1):end])
                    tentative_index -= 1
                    reduc = 1
                end

                if NEAR_SINGULAR[] || dt * T_unit < 1e-10
                    if reduc == 1
                        tentative_index += 1
                    end
                    tvec = filter(x -> !(x >= tvec[tentative_index] && x < (tvec[tentative_index] + dtb)), tvec)
                    g = (a, b) -> a + 0.6 * (b - a)
                    tvec = vcat(tvec[1:(tentative_index - 1)],
                        g(tvec[tentative_index - 1], tvec[tentative_index]),
                        tvec[tentative_index:end])
                    tentative_index -= 2
                    reduc = 1

                    ps_accepted = ps_old
                    previous_conditions[2] = previous_conditions[1]
                    previous_conditions[1] = old_conditions
                    going_back += 1
                    println("Warning detected. Will proceed with going back")
                    if going_back > 50
                        error("Went back too many times. Stopping execution")
                    end
                elseif reduc == 0
                    if norm(psNew, 1) == 0
                        B_l_ps_new = zeros(Float64, N)
                    else
                        nb_contact_points = Int(round(numlTent))
                        if PROBLEM_CONSTANTS["linear_on_theta"]
                            if nb_contact_points > 1
                                contactAngle = 1.5 * thetaVec[nb_contact_points] - 0.5 * thetaVec[nb_contact_points - 1]
                            else
                                contactAngle = (thetaVec[2] + thetaVec[1]) / 2
                            end
                            angles = range(contactAngle, stop=pi, length=(nb_contact_points + 1) * PROBLEM_CONSTANTS["interpolation_number"])
                            values = r_from_spherical(collect(angles), oscillation_amplitudes[:, tentative_index])
                            f_interp = r -> interp1_linear(dr .* collect(0:nb_contact_points),
                                vcat(psNew[1:nb_contact_points], 0.0), r, extrap=0.0)
                            values = f_interp.(values)
                            B_l_ps_new = custom_project_amplitudes(collect(angles), values, N, NaN, NaN)
                        else
                            error("project_amplitudes path is not implemented in this translation.")
                        end
                    end

                    amplitudes_new, velocities_new = solve_ODE_unkown(NaN, B_l_ps_new, dt, previous_conditions, PROBLEM_CONSTANTS)
                    errorP = norm(amplitudes_tent - amplitudes_new) / norm(amplitudes_tent)

                    if PROBLEM_CONSTANTS["DEBUG_FLAG"]
                        println(@sprintf("Inside ll: %0.2g, errP: %0.5g", ll, errorP))
                    end

                    if errorP < tolP
                        ps_old = ps_accepted
                        old_conditions = previous_conditions[1]

                        numl[tentative_index + 1] = numlTent
                        eta_accepted = etaTent
                        phi_accepted = phiTent
                        ps_accepted = psNew
                        z[tentative_index + 1] = zTent
                        vz[tentative_index + 1] = vzTent
                        oscillation_amplitudes[:, tentative_index + 1] .= amplitudes_new
                        pressure_amplitudes[:, tentative_index + 1] .= B_l_ps_new
                        Rv[tentative_index + 1] = zs_from_spherical(pi, amplitudes_new)

                        previous_conditions[1] = previous_conditions[2]
                        previous_conditions[2] = Condition(copy(amplitudes_new), copy(velocities_new), copy(B_l_ps_new), dt, N,
                            previous_conditions[1].current_time + dt, z[tentative_index + 1], vz[tentative_index + 1], numlTent)

                        nlmax[tentative_index + 1] = nlmaxTent
                        etaOri[tentative_index + 1] = eta_accepted[1]
                        etas[:, tentative_index + 1] .= eta_accepted
                        phis[:, tentative_index + 1] .= phi_accepted

                        if current_to_save <= length(tvecOri) && (t >= tvecOri[current_to_save] || t < 10 * dtb)
                            indexes_to_save[current_to_save] = tentative_index
                            current_to_save += 1
                        end

                        if (zTent > 1.5 && numlTent == 0) || (vz[tentative_index + 1] < 0 && vz[tentative_index] >= 0 && t > 50 * dtb)
                            tend = t
                        end
                    else
                        amplitudes_tent = amplitudes_new
                        RmaxTent = r_from_spherical(maximum_contact_radius(amplitudes_tent), amplitudes_tent)
                        nlmaxTent = floor(RmaxTent / dr) + 1
                        nlmaxTent_int = Int(round(nlmaxTent))
                        thetaVec = theta_from_cylindrical(dr .* collect(0:(nlmaxTent_int - 1)),
                            oscillation_amplitudes[:, tentative_index])
                        RvTent = zs_from_spherical(pi, amplitudes_tent)
                        zs[1:nlmaxTent_int] .= zs_from_spherical(thetaVec, amplitudes_tent) .- RvTent
                        zs[(nlmaxTent_int + 1):nr] .= Inf

                        tanDrop = calculate_tan(dr .* collect(1:nlmaxTent_int) .- dr / 2, amplitudes_tent)
                        angleDropMP[1:nlmaxTent_int] .= atan.(tanDrop[1:nlmaxTent_int])
                        psprob = zeros(Float64, nlmaxTent_int, 5)
                    end
                end
            end
        end

        max_idx = min(current_to_save - 1, length(indexes_to_save))
        indexes_to_save = normalize_indexes(indexes_to_save[1:max_idx])
        results = Dict(
            "z" => z[indexes_to_save],
            "etaOri" => etaOri[indexes_to_save],
            "etas" => etas[:, indexes_to_save],
            "vz" => vz[indexes_to_save],
            "tvec" => tvec[indexes_to_save],
            "nlmax" => nlmax[indexes_to_save],
            "numl" => numl[indexes_to_save],
            "oscillation_amplitudes" => oscillation_amplitudes[:, indexes_to_save],
            "pressure_amplitudes" => pressure_amplitudes[:, indexes_to_save],
            "Rv" => Rv[indexes_to_save],
        )
        save_results(joinpath(out_dir, prefix * "results.jld2"), results)
    catch err
        max_idx = min(current_to_save - 1, length(indexes_to_save))
        indexes_to_save = normalize_indexes(indexes_to_save[1:max_idx])
        results = Dict(
            "z" => z[indexes_to_save],
            "etaOri" => etaOri[indexes_to_save],
            "etas" => etas[:, indexes_to_save],
            "vz" => vz[indexes_to_save],
            "tvec" => tvec[indexes_to_save],
            "nlmax" => nlmax[indexes_to_save],
            "numl" => numl[indexes_to_save],
            "oscillation_amplitudes" => oscillation_amplitudes[:, indexes_to_save],
            "pressure_amplitudes" => pressure_amplitudes[:, indexes_to_save],
            "Rv" => Rv[indexes_to_save],
        )
        save_results(joinpath(out_dir, prefix * "errored_results.jld2"), results)

        println(@sprintf("Couldn't run simulation with the following parameters:\n Velocity: %g \n Modes: %g", U0, N))
        stamp = Dates.format(now(), "yyyymmddMMSS")
        open(joinpath(out_dir, @sprintf("error_logU0=%g-%s.txt", U0, stamp)), "w") do io
            write(io, string(err))
        end
        rethrow(err)
    end

    simul_time = time() - tstart
    conditions = Dict(
        "T" => T_unit,
        "N" => N,
        "U0" => U0,
        "Ang" => Ang,
        "Re" => Re,
        "Fr" => Fr,
        "We" => We,
        "WeS" => WeS,
        "Cang" => Cang,
        "tend" => tend,
        "nsteps" => nsteps,
        "dtb" => dtb,
        "L_unit" => L_unit,
        "T_unit" => T_unit,
        "M_unit" => M_unit,
        "PROBLEM_CONSTANTS" => PROBLEM_CONSTANTS,
        "simul_time" => simul_time,
    )
    save_results(joinpath(out_dir, prefix * "problem_conditions.jld2"), conditions)
    println(@sprintf("Finished simulation. Time elapsed: %0.2f minutes", simul_time / 60))
end
