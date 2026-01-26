function solveDD0(dt, z, vz, etao, phio, nr, Re, Delta, DTN, Fr, We, zs, Rv)
    b = vcat(etao, phio)
    Sist = [eye(nr) - dt * 2 * Delta / Re  -dt * DTN;
            dt * (eye(nr) / Fr - Delta / We)  eye(nr) - dt * 2 * Delta / Re]
    c = solve_linear_system(Sist, b)
    etaprob0 = c[1:nr]
    phiprob0 = c[(nr + 1):(2 * nr)]
    vzprob0 = vz - dt / Fr
    zprob0 = z + vz * dt - dt^2 / (2 * Fr)
    check = etaprob0 .> (zprob0 .+ zs .+ Rv)
    errortan = sum(check) > 0.5 ? 4.0 : 0.0
    return etaprob0, phiprob0, zprob0, vzprob0, errortan
end

function solvenDDCusp(nlprev, nl, dt, zo, vzo, etao, phio, nr, dr, Re, Delta, DTN, Fr, We, _Ma, zs, Int, angleDrop, Cang, Dr, Rv)
    Ma = 4 * pi / 3
    We = We / Dr
    b = vcat(etao, phio)
    if nlprev < nl
        AngC = min(angleDrop[nl], pi - Cang)
    else
        AngC = pi - Cang
    end

    Sist = [eye(nr) - dt * 2 * Delta / Re  -dt * DTN;
            dt * (eye(nr) / Fr - Delta / We)  eye(nr) - dt * 2 * Delta / Re]
    bmod = b - Sist[:, 1:nl] * (zs[1:nl] .+ Rv)

    block_ps = vcat(zeros(nr, nl), dt * eye(nl) * Dr, zeros(nr - nl, nl))
    int_row = reshape(Int[1:nl], 1, :)
    Mat = [Sist[:, (nl + 1):(2 * nr)]  block_ps  zeros(2 * nr, 1)  Sist[:, 1:nl] * ones(nl);
           zeros(1, 2 * nr - nl)  -dt * int_row / Ma  1  0;
           zeros(1, 2 * nr - nl)  zeros(1, nl)  -dt  1]

    indep = vcat(bmod, vzo - dt / Fr, zo)
    ds = solve_linear_system(Mat, indep)

    etaprob = zeros(Float64, nr)
    etaprob[(nl + 1):nr] .= ds[1:(nr - nl)]
    phiprob = ds[(nr - nl + 1):(2 * nr - nl)]
    psprob = ds[(2 * nr - nl + 1):(2 * nr)]
    vzprob = ds[end - 1]
    zprob = ds[end]
    etaprob[1:nl] .= zprob .+ Rv .+ zs[1:nl]

    check = etaprob[(nl + 1):end] .> (zprob .+ Rv .+ zs[(nl + 1):end])
    if sum(check) > 0.5
        errortan = 4.0
    else
        taneffect = (etaprob[nl + 1] - etaprob[nl]) / dr
        ataneffect = atan(taneffect)
        errortan = angleDrop[nl] - ataneffect - AngC
    end
    return etaprob, phiprob, zprob, vzprob, psprob, errortan
end

