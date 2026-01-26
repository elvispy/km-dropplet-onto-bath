function collectPl(n::Int, x)
    xv = vec(x)
    y = ones(Float64, n + 1, length(xv))
    if n >= 1
        y[2, :] .= xv
    end
    for idx in 3:(n + 1)
        y[idx, :] .= ((2 * idx - 3) / (idx - 1)) .* xv .* y[idx - 1, :] .-
                     ((idx - 2) / (idx - 1)) .* y[idx - 2, :]
    end
    return y[2:end, :]
end

function collectdnPl(nmax::Int, x, order::Int=1)
    xv = vec(x)
    if order == 1
        y = ones(Float64, nmax, length(xv))
        if nmax >= 2
            y[2, :] .= 3 .* xv
        end
        for idx in 2:(nmax - 1)
            y[idx + 1, :] .= ((2 * idx + 1) / idx) .* xv .* y[idx, :] .-
                             ((idx + 1) / idx) .* y[idx - 1, :]
        end
    else
        y = zeros(Float64, nmax, length(xv))
        if nmax >= 2
            y[2, :] .= 3
        end
        if nmax >= 3
            y[3, :] .= 15 .* xv
        end
        for idx in 3:(nmax - 1)
            y[idx + 1, :] .= ((2 * idx + 1) / (idx - 1)) .* xv .* y[idx, :] .-
                             ((idx + 2) / (idx - 1)) .* y[idx - 1, :]
        end
    end
    return y
end

function legendre_sum(amplitudes::AbstractVector{<:Real}, cosθ)
    cos_vec = cosθ isa Number ? [cosθ] : vec(cosθ)
    pl = collectPl(length(amplitudes), cos_vec)
    vals = vec(sum(amplitudes .* pl, dims=1))
    return cosθ isa Number ? vals[1] : reshape(vals, size(cosθ))
end

function legendre_dn_sum(amplitudes::AbstractVector{<:Real}, cosθ, order::Int=1)
    cos_vec = cosθ isa Number ? [cosθ] : vec(cosθ)
    pl = collectdnPl(length(amplitudes), cos_vec, order)
    vals = vec(sum(amplitudes .* pl, dims=1))
    return cosθ isa Number ? vals[1] : reshape(vals, size(cosθ))
end

function zeta_generator(amplitudes)
    amps = amplitudes isa AbstractArray ? vec(amplitudes) : [amplitudes]
    return theta -> legendre_sum(amps, cos.(theta))
end
