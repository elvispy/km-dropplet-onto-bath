function build_domain(D::Real, quant::Int)
    nr = ceil(Int, D * quant / 2)
    dr = D / (2 * nr)
    nlmax = quant + 1

    r = collect(0:dr:(D / 2))
    rn = collect(0:(nr + 1))
    xplot = vcat(-reverse(r[2:(nr + 1)]), r)

    zs_vals = zeros(Float64, nlmax)
    zs_vals[1] = 0.0
    @inbounds for j in 1:quant
        zs_vals[j + 1] = 1 - sqrt(1 - (j * dr)^2)
    end
    zs = if nr > nlmax
        vcat(zs_vals, fill(100.0, nr - nlmax))
    else
        zs_vals[1:nr]
    end

    subdiag = zeros(Float64, nr, nr)
    @inbounds for i in 2:nr
        subdiag[i, i - 1] = 1.0
    end

    Dr = (-subdiag + subdiag') / 2
    OneOverR = zeros(Float64, nr)
    @inbounds for i in 2:nr
        OneOverR[i] = 1 / r[i]
    end

    DrOverR = zeros(Float64, nr, nr)
    @inbounds for i in 2:nr
        DrOverR[i, :] .= OneOverR[i] .* Dr[i, :]
    end

    Derivr = Dr / dr
    Derivr[1, :] .= 0.0

    Derivrr = (subdiag + subdiag' - 2 .* eye(nr)) / dr^2
    Derivrr[1, :] .= 0.0
    Derivrr[1, 1] = -2 / dr^2
    if nr >= 2
        Derivrr[1, 2] = 2 / dr^2
    end

    Delta = Derivrr + DrOverR / dr
    Delta[1, :] .= 0.0
    Delta[1, 1] = -4 / dr^2
    if nr >= 2
        Delta[1, 2] = 4 / dr^2
    end

    tanDropMP = zeros(Float64, nlmax)
    @inbounds for kk in 1:(nlmax - 1)
        tanDropMP[kk] = dr * (kk - 0.5) / sqrt(1 - (kk - 0.5)^2 * dr^2)
    end
    angleDropMP = atan.(tanDropMP[1:(nlmax - 1)])
    angleDropMP = vcat(angleDropMP, pi / 2)

    nint = 2 * nlmax
    IntMat = zeros(Float64, nint, nint)
    IntMat[1, 1] = 1 / 12
    @inbounds for ii in 2:nint
        IntMat[ii, 1] = 1 / 3
        for jj in 2:(ii - 1)
            IntMat[ii, jj] = 2 * jj - 2
        end
        IntMat[ii, ii] = 1.5 * ii - 21 / 12
    end
    IntMat .*= pi * dr^2

    return (
        D = float(D),
        quant = quant,
        nr = nr,
        dr = dr,
        Delta = Delta,
        IntMat = IntMat,
        rn = rn,
        r = r,
        xplot = xplot,
        zs = zs,
        angleDropMP = angleDropMP,
    )
end
