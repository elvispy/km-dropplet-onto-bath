function r_from_spherical(angles, amplitudes)
    amps = amplitudes isa AbstractArray ? vec(amplitudes) : [amplitudes]
    zeta = legendre_sum(amps, cos.(angles))
    return sin.(angles) .* (1 .+ zeta)
end

function zs_from_spherical(theta, amplitudes)
    amps = amplitudes isa AbstractArray ? vec(amplitudes) : [amplitudes]
    zeta = legendre_sum(amps, cos.(theta))
    return cos.(theta) .* (1 .+ zeta)
end

function theta_from_cylindrical(r, amplitudes)
    amps = amplitudes isa AbstractArray ? vec(amplitudes) : [amplitudes]
    zeta = zeta_generator(amps)
    f_prime = theta -> cos(theta) .* (1 .+ zeta(theta)) .-
        sin(theta).^2 .* legendre_dn_sum(amps, cos(theta), 1)

    angle = zeros(Float64, length(r))
    rvec = vec(r)
    for ii in eachindex(rvec)
        if rvec[ii] == 0
            angle[ii] = pi
            continue
        end
        f_objective = theta -> sin(theta) .* (1 .+ zeta(theta)) .- rvec[ii]

        if ii >= 5
            theta = interp1_makima(rvec[(ii - 4):(ii - 1)], angle[(ii - 4):(ii - 1)], rvec[ii])
        else
            theta = pi - asin(min(1, rvec[ii]))
        end

        tol_theta = 1e-7
        n = 1
        while abs(f_objective(theta)) >= tol_theta && n < 350
            theta = mod(theta - f_objective(theta) / f_prime(theta) - 1e-4, pi / 2) + 1e-4 + pi / 2
            n += 1
            if n == 200
                theta = pi
            elseif n == 300
                theta = pi - asin(min(1, rvec[ii]))
            end
        end
        if n >= 340
            @warn "Theta_from_cylindrical didnt converge"
        end
        angle[ii] = theta
    end
    return angle
end

function calculate_tan(distance_to_axis, amplitudes)
    amps = amplitudes isa AbstractArray ? vec(amplitudes) : [amplitudes]
    angle = theta_from_cylindrical(distance_to_axis, amps)
    zeta = zeta_generator(amps)
    der = theta -> legendre_dn_sum(amps, cos(theta), 1)
    dzdr = theta -> (-sin(theta) .* (1 .+ zeta(theta)) .-
                     cos(theta) .* sin(theta) .* der(theta)) ./
        (cos(theta) .* (1 .+ zeta(theta)) .- sin(theta).^2 .* der(theta))
    return tan.(dzdr.(angle))
end

function maximum_contact_radius(amplitudes)
    amps = amplitudes isa AbstractArray ? vec(amplitudes) : [amplitudes]
    zeta = zeta_generator(amps)
    drdtheta = theta -> cos(theta) * (1 + zeta(theta)) - sin(theta)^2 * legendre_dn_sum(amps, cos(theta), 1)
    dr2dtheta2 = theta -> -sin(theta) * (1 + zeta(theta)) - 2 * cos(theta) * sin(theta) *
        legendre_dn_sum(amps, cos(theta), 1) + sin(theta)^3 * legendre_dn_sum(amps, cos(theta), 2)

    theta = pi / 2 + pi / 4
    tol_theta = 1e-7
    n = 1
    while abs(drdtheta(theta)) >= tol_theta && n < 150
        theta = mod(theta - drdtheta(theta) / dr2dtheta2(theta) - 1e-4, pi) + 1e-4
        n += 1
        if n == 50
            theta = pi / 2
        elseif n == 100
            theta = rand() / 100 + pi / 2
        end
        if n == 149
            thetas = range(pi / 3, 4 * pi / 5, length=1000)
            return maximum(r_from_spherical(collect(thetas), amps))
        end
    end
    return r_from_spherical(theta, amps)
end
