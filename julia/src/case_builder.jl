function build_case(; D=5.0, quant=20, rho=1.0, sigma=72.2, nu=9.78e-3, muair=0.0, g=981.0,
    rhoS=1.0, sigmaS=72.2, Ro=0.035, refp=10, num_batches=20, threaded=true, output_path=nothing)

    domain = build_domain(D, quant)
    DTN = generate_dtn_new345(domain.nr, D; refp=refp, num_batches=num_batches, threaded=threaded)

    Ma = 4 * pi * rhoS / (3 * rho)
    Ra = rhoS / rho

    case = Dict(
        "Ro" => float(Ro),
        "rhoS" => float(rhoS),
        "sigmaS" => float(sigmaS),
        "rho" => float(rho),
        "sigma" => float(sigma),
        "nu" => float(nu),
        "muair" => float(muair),
        "g" => float(g),
        "D" => float(D),
        "quant" => Int(quant),
        "nr" => Int(domain.nr),
        "dr" => float(domain.dr),
        "Delta" => domain.Delta,
        "IntMat" => domain.IntMat,
        "DTN" => DTN,
        "Ma" => float(Ma),
        "Ra" => float(Ra),
    )

    if output_path !== nothing
        save_results(string(output_path), case)
    end

    return case
end
