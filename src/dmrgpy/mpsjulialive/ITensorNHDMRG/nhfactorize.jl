function nhfactorize(
    ::Algorithm"biorthoblock", phil, phir, drho, lindsl, lindsr, targettags, identities; kwargs...
)
    # https://arxiv.org/abs/2401.15000
    # compute reduced density matrix and apply perturbation
    rho = dag(replaceinds(phil, lindsl, lindsl')) 
    if !isnothing(identities)
        rho *= identities[3]
        
        cinds = commoninds(identities[3], phil)
        rho = replaceinds(rho, cinds', cinds)
    end
    rho *= phir

    if !isnothing(drho)
        rho += swapprime(drho, 0=>1)
    end

    B, Y, Ybar, spec = biorthoblocktransform(rho, lindsl', dag(lindsr); kwargs...)
    noprime!(Y)
    noprime!(Ybar)
    Y = replaceinds!(Y, commonind(Y, B) => settags(commonind(Y, B), targettags))
    Ybar = replaceinds!(Ybar, commonind(Ybar, B) => settags(commonind(Ybar, B), targettags))

    return dag(Y), Ybar, spec
end

function nhfactorize(
    ::Algorithm"fidelity", phil, phir, drho, lindsl, lindsr, targettags, identities; kwargs...
)
    # Phys. Rev. B 105, 205125 
    # https://doi.org/10.1103/PhysRevB.105.205125
    # compute reduced density matrix and apply perturbation
    idleft, idright = nothing, nothing 
    if !isnothing(identities)
        idleft, idright, _ = identities
    end

    rhor = replaceinds(phir, lindsr, lindsr')
    if !isnothing(idright)
        rhor = prime(rhor, commonind(rhor, idright))
        rhor *= idright
    end
    rhor *= dag(phir)
    
    rhol = replaceinds(phil, lindsl, lindsl')
    if !isnothing(idleft)
        rhol = prime(rhol, commonind(rhol, idleft))
        rhol *= idleft
    end
    rhol *= dag(phil)
 
    !hassameinds(rhol, rhor) && error("Left and right states need to share the same indices")

    rho = (rhol + rhor) / 2
    if !isnothing(drho)
        rho += drho
    end

    D, U, spec = eigen(rho, lindsl', dag(lindsl); ishermitian=true, kwargs...)
    U = noprime!(U)
    return U, U, spec
end

function nhreplacebond!(
    Ml::MPS,
    Mr::MPS,
    b::Int,
    phil::ITensor,
    phir::ITensor,
    alg,
    identities=nothing;
    ortho=nothing,
    eigen_perturbation=nothing,
    # Decomposition kwargs
    mindim=nothing,
    maxdim=nothing,
    cutoff=nothing,
    kwargs...,
)
    ortho = NDTensors.replace_nothing(ortho, "left")

    if ortho != "left" && ortho != "right"
        error(
            "In replacebond!, got ortho = $ortho, only currently supports `left` and `right`.",
        )
    end

    leftindsl = if ortho == "left"
        commoninds(Ml[b], phil)
    else
        commoninds(Ml[b+1], phil)
    end
    leftindsr = if ortho == "left"
        commoninds(Mr[b], phir)
    else
        commoninds(Mr[b+1], phir)
    end

    U, Ubar, spec = nhfactorize(
        Algorithm(alg),
        phil,
        phir,
        eigen_perturbation,
        leftindsl,
        leftindsr,
        tags(commonind(Ml[b], Ml[b+1])),
        identities;
        mindim,
        maxdim,
        cutoff,
        kwargs...,
    )

    for (M, phi, U1, U2) in [(Ml, phil, Ubar, U), (Mr, phir, U, Ubar)]
        L, R = if ortho == "left"
            U1, phi * dag(U2)
        elseif ortho == "right"
            phi * dag(U2), U1
        end
        M[b] = L
        M[b+1] = R
    end

    return spec
end

function biorthogonalize!(psil, psir, alg; mindim=nothing, maxdim=10, cutoff=nothing, kwargs...)
    @assert siteinds(psir) == siteinds(psil) "both MPS need to share the same basis"

    @assert abs(inner(psil, psir)) >= sqrt(cutoff) "The initial vectors are almost orthogonal, overlap is $(inner(psil, psir))"

    sites = siteinds(psir)

    noprime!(psir)
    noprime!(psil)

    function buildidentities(psi1, psi2)
        M = ITensor[]
        for i in firstindex(sites):lastindex(sites)
            Mi = length(M) > 0 ? last(M) : ITensor(1)

            Mi *= prime(psi1[i], "Link")
            Mi *= dag(psi2[i])
            push!(M, Mi)
        end
        return M
    end

    Mll = buildidentities(psil, psil)
    Mrr = buildidentities(psir, psir)
    Mlr = buildidentities(psir, psil)

    for i in (lastindex(sites)-1):-1:1
        phil = psil[i] * psil[i+1]
        phir = psir[i] * psir[i+1]

        identities = nothing 
        if i > 1
            identities = (Mll[i-1], Mrr[i-1], Mlr[i-1])
        end

        nhreplacebond!(
            psil, psir, i, phil, phir, alg, identities; ortho="right", mindim, maxdim, cutoff, kwargs...
        )
    end

    noprime!(psil)
    noprime!(psir)

    return psil, psir
end
