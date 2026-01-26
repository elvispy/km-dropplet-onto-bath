function custom_project_amplitudes(angles, values, N, _unused=nothing, _unused2=nothing)
    cosang = cos.(angles)
    legendre_values = vcat(ones(1, length(cosang)), collectPl(N + 2, cosang))
    M = collect(3:2:(2 * N + 3))
    M = reshape(M, :, 1)

    Gn = (legendre_values[3:end, :] .- legendre_values[1:(end - 2), :]) ./ M
    Gn = vcat(reshape(cosang, 1, :), Gn)
    Gnab = diff(Gn, dims=2)

    idx_hi = reshape(collect(2:(N + 1)), :, 1)
    idx_lo = reshape(collect(1:N), :, 1)
    intxPn = (idx_hi .* Gnab[3:end, :] .+ idx_lo .* Gnab[1:N, :]) ./ M[1:N]

    as = diff(values) ./ diff(cosang)
    bs = values[1:(end - 1)] .- as .* cosang[1:(end - 1)]
    as_row = reshape(as, 1, :)
    bs_row = reshape(bs, 1, :)

    coefs = -M[1:(end - 1)] ./ 2 .* sum(as_row .* intxPn .+ bs_row .* Gnab[2:(end - 1), :], dims=2)
    return vec(coefs)
end

