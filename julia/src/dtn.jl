function _batch_ranges(total::Int, num_batches::Int)
    if total <= 0
        return UnitRange{Int}[]
    end
    batches = UnitRange{Int}[]
    if num_batches < 1
        push!(batches, 1:total)
        return batches
    end
    batch_size = ceil(Int, total / num_batches)
    start_idx = 1
    while start_idx <= total
        stop_idx = min(start_idx + batch_size - 1, total)
        push!(batches, start_idx:stop_idx)
        start_idx = stop_idx + 1
    end
    return batches
end

function generate_dtn_new345(nr::Int, D::Real; refp::Int=10, num_batches::Int=20, threaded::Bool=true)
    dr = D / (2 * nr)
    rn = collect(0:(nr + 1))
    drp = dr / refp

    numer = ceil(Int, pi * D / drp)
    if isodd(numer)
        numer += 1
    end
    dtheta = 2 * pi / numer

    l_vals = collect(range(dtheta / 2, step=dtheta, stop=pi - dtheta / 4))
    cos_l = cos.(l_vals)
    sin_l = sin.(l_vals)
    batches = _batch_ranges(length(l_vals), num_batches)

    DTN = zeros(Float64, nr, nr)

    # k = 1: integrate away from singularity using analytic weights
    row1 = zeros(Float64, nr)
    rn_k = Float64(rn[1])
    i_end = Int(round(rn_k + nr + 1))
    @inbounds for i in 2:i_end
        idx1 = Int(round(i))
        term0 = (-i^2 / 2 - i - 1 / 3) * log((i + 1) / i) +
            (i^2 / 6 + i / 2 + 1 / 3) * (1 - i / (i + 1)) -
            (i + 0.5) / 6 + i / 2 + 1 / 2
        term1 = (3 * i^2 / 2 + 2 * i - 1 / 2) * log((i + 1) / i) +
            (-i^2 / 2 - i + 1 / 2 + 1 / i) * (1 - i / (i + 1)) +
            (i + 0.5) / 2 - 3 * i / 2 - 1
        term2 = (-3 * i^2 / 2 - i + 1) * log((i + 1) / i) +
            (i^2 / 2 + i / 2 - 1) * (1 - i / (i + 1)) -
            (i + 0.5) / 2 + 3 * i / 2 + 1 / 2
        term3 = (i^2 / 2 - 1 / 6) * log((i + 1) / i) +
            (-i^2 + 1) / 6 * (1 - i / (i + 1)) +
            (i + 0.5) / 6 - i / 2

        if idx1 <= nr
            row1[idx1] -= term0
        end
        if idx1 < nr
            row1[idx1 + 1] -= term1
            if idx1 < nr - 1
                row1[idx1 + 2] -= term2
                if idx1 < nr - 2
                    row1[idx1 + 3] -= term3
                end
            end
        end
    end

    row1[1] += 1 / 2
    row1 ./= dr
    row1[1] = row1[1] + 209 / (54 * dr)
    if nr >= 2
        row1[2] = row1[2] - 29 / (6 * dr)
    end
    if nr >= 3
        row1[3] = row1[3] + 7 / (6 * dr)
    end
    if nr >= 4
        row1[4] = row1[4] - 11 / (54 * dr)
    end
    DTN[1, :] .= row1

    function integrate_far!(line::Vector{Float64}, rn_k::Float64, cosb::AbstractVector{Float64},
        sinb::AbstractVector{Float64}, i_start::Int, i_end::Int)
        if i_end < i_start
            return
        end
        @inbounds for i in i_start:i_end
            kern = 2 * (1 / (i - 0.5) - 1 / (i + 0.5))
            for jj in eachindex(cosb)
                cosl = cosb[jj]
                sinl = sinb[jj]
                radn = sqrt((rn_k + i * cosl / refp)^2 + (i * sinl / refp)^2)
                idx = floor(Int, radn)
                w1 = radn - idx
                if w1 < 0
                    w1 = 0.0
                elseif w1 > 1
                    w1 = 1.0
                end
                w2 = w1 * w1
                w3 = w2 * w1

                if idx == 0
                    line[1] -= (3 * w3 / 4 - 7 * w2 / 4 + 1) * kern
                    if nr >= 2
                        line[2] -= (-w3 + 2 * w2) * kern
                    end
                    if nr >= 3
                        line[3] -= (w3 - w2) / 4 * kern
                    end
                elseif idx < nr
                    line[idx] -= (-w3 / 6 + w2 / 2 - w1 / 3) * kern
                    if idx + 1 <= nr
                        line[idx + 1] -= (w3 / 2 - w2 - w1 / 2 + 1) * kern
                    end
                    if idx + 2 <= nr - 0
                        if idx < nr - 1
                            line[idx + 2] -= (-w3 / 2 + w2 / 2 + w1) * kern
                        end
                        if idx < nr - 2
                            line[idx + 3] -= (w3 - w1) / 6 * kern
                        end
                    end
                end
            end
        end
    end

    # k = 2 and k = 3 (special near-origin coefficients)
    for k in 2:min(3, nr)
        line = zeros(Float64, nr)
        rn_k = Float64(rn[k])
        for batch in batches
            cosb = @view cos_l[batch]
            sinb = @view sin_l[batch]

            @inbounds for i in 1:(2 * refp)
                kern = 2 * (1 / (i - 0.5) - 1 / (i + 0.5))
                for jj in eachindex(cosb)
                    cosl = cosb[jj]
                    sinl = sinb[jj]
                    radn = sqrt((rn_k + i * cosl / refp)^2 + (i * sinl / refp)^2)
                    x1 = i * cosl / refp
                    posr = radn - rn_k
                    posr2 = posr * posr
                    posr3 = posr2 * posr
                    posr4 = posr3 * posr

                    if k == 2
                        line[k] -= ((-6 - 3 * posr + 9 * posr2 + 3 * posr3 - 3 * posr4) * posr + 6 * x1) / 72 * kern
                        line[k - 1] -= ((-88 + 48 * posr + 62 * posr2 - 12 * posr3 - 10 * posr4) * posr + 88 * x1) / 72 * kern
                        line[k] -= ((72 - 90 * posr - 90 * posr2 + 18 * posr3 + 18 * posr4) * posr - 72 * x1) / 72 * kern
                        if k <= nr - 1
                            line[k + 1] -= ((24 + 48 * posr + 18 * posr2 - 12 * posr3 - 6 * posr4) * posr - 24 * x1) / 72 * kern
                        end
                        if k <= nr - 2
                            line[k + 2] -= ((-2 - 3 * posr + posr2 + 3 * posr3 + posr4) * posr + 2 * x1) / 72 * kern
                        end
                    else
                        line[k - 2] -= ((124 - 12 * posr - 149 * posr2 + 12 * posr3 + 25 * posr4) * posr - 124 * x1) / 288 * kern
                        line[k - 1] -= ((-384 + 192 * posr + 288 * posr2 - 48 * posr3 - 48 * posr4) * posr + 384 * x1) / 288 * kern
                        line[k] += ((-4 + 10 * posr + 5 * posr2 - 2 * posr3 - posr4) * posr + 4 * x1) / 8 * kern
                        if k <= nr - 1
                            line[k + 1] -= ((128 + 192 * posr + 32 * posr2 - 48 * posr3 - 16 * posr4) * posr - 128 * x1) / 288 * kern
                        end
                        if k <= nr - 2
                            line[k + 2] -= ((-12 - 12 * posr + 9 * posr2 + 12 * posr3 + 3 * posr4) * posr + 12 * x1) / 288 * kern
                        end
                    end
                end
            end

            i_start = 2 * refp + 1
            i_end = Int(round((rn_k + nr) * refp))
            integrate_far!(line, rn_k, cosb, sinb, i_start, i_end)
        end

        line .*= dtheta / (2 * pi * drp)
        line[k] += 2 / (4 * dr + drp)
        DTN[k, :] .= line
    end

    # k >= 4
    if nr >= 4
        if threaded && Threads.nthreads() > 1
            Threads.@threads for k in 4:nr
                line = zeros(Float64, nr)
                rn_k = Float64(rn[k])
                for batch in batches
                    cosb = @view cos_l[batch]
                    sinb = @view sin_l[batch]

                    @inbounds for i in 1:(2 * refp)
                        kern = 2 * (1 / (i - 0.5) - 1 / (i + 0.5))
                        for jj in eachindex(cosb)
                            cosl = cosb[jj]
                            sinl = sinb[jj]
                            radn = sqrt((rn_k + i * cosl / refp)^2 + (i * sinl / refp)^2)
                            x1 = i * cosl / refp
                            posr = radn - rn_k
                            posr2 = posr * posr
                            posr3 = posr2 * posr

                            line[k - 2] -= ((2 - posr - 2 * posr2 + posr3) * posr - 2 * x1) / 24 * kern
                            line[k - 1] -= (4 * (-4 + 4 * posr + posr2 - posr3) * posr + 16 * x1) / 24 * kern
                            line[k] -= (posr2 - 5) * posr2 / 4 * kern
                            if k <= nr - 1
                                line[k + 1] -= (4 * (4 + 4 * posr - posr2 - posr3) * posr - 16 * x1) / 24 * kern
                            end
                            if k <= nr - 2
                                line[k + 2] -= ((-2 - posr + 2 * posr2 + posr3) * posr + 2 * x1) / 24 * kern
                            end
                        end
                    end

                    i_start = 2 * refp + 1
                    i_end = Int(round((rn_k + nr) * refp))
                    integrate_far!(line, rn_k, cosb, sinb, i_start, i_end)
                end

                line .*= dtheta / (2 * pi * drp)
                line[k] += 2 / (4 * dr + drp)
                DTN[k, :] .= line
            end
        else
            for k in 4:nr
                line = zeros(Float64, nr)
                rn_k = Float64(rn[k])
                for batch in batches
                    cosb = @view cos_l[batch]
                    sinb = @view sin_l[batch]

                    @inbounds for i in 1:(2 * refp)
                        kern = 2 * (1 / (i - 0.5) - 1 / (i + 0.5))
                        for jj in eachindex(cosb)
                            cosl = cosb[jj]
                            sinl = sinb[jj]
                            radn = sqrt((rn_k + i * cosl / refp)^2 + (i * sinl / refp)^2)
                            x1 = i * cosl / refp
                            posr = radn - rn_k
                            posr2 = posr * posr
                            posr3 = posr2 * posr

                            line[k - 2] -= ((2 - posr - 2 * posr2 + posr3) * posr - 2 * x1) / 24 * kern
                            line[k - 1] -= (4 * (-4 + 4 * posr + posr2 - posr3) * posr + 16 * x1) / 24 * kern
                            line[k] -= (posr2 - 5) * posr2 / 4 * kern
                            if k <= nr - 1
                                line[k + 1] -= (4 * (4 + 4 * posr - posr2 - posr3) * posr - 16 * x1) / 24 * kern
                            end
                            if k <= nr - 2
                                line[k + 2] -= ((-2 - posr + 2 * posr2 + posr3) * posr + 2 * x1) / 24 * kern
                            end
                        end
                    end

                    i_start = 2 * refp + 1
                    i_end = Int(round((rn_k + nr) * refp))
                    integrate_far!(line, rn_k, cosb, sinb, i_start, i_end)
                end

                line .*= dtheta / (2 * pi * drp)
                line[k] += 2 / (4 * dr + drp)
                DTN[k, :] .= line
            end
        end
    end

    return DTN
end
