function interp1_linear(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, xq; extrap=0.0)
    if xq isa AbstractVector
        return [interp1_linear(x, y, xi; extrap=extrap) for xi in xq]
    end
    if xq < x[1] || xq > x[end]
        return extrap
    end
    idx = searchsortedlast(x, xq)
    if idx >= length(x)
        return y[end]
    end
    x0 = x[idx]
    x1 = x[idx + 1]
    y0 = y[idx]
    y1 = y[idx + 1]
    if x1 == x0
        return y0
    end
    return y0 + (y1 - y0) * ((xq - x0) / (x1 - x0))
end

function interp1_makima(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, xq)
    return interp1_linear(x, y, xq; extrap=y[end])
end
