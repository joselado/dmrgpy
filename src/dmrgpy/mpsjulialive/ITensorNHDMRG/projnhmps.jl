using ITensors
import ITensorMPS: MPS

mutable struct ProjNHMPS
    projL::ITensorMPS.ProjMPS
    projR::ITensorMPS.ProjMPS
    cached_site_range::AbstractRange
    envL::ITensor
    envR::ITensor
end
function ProjNHMPS(Ml::MPS, Mr::MPS)
    return ProjNHMPS(
        ITensorMPS.ProjMPS(Ml),
        ITensorMPS.ProjMPS(Mr),
        -1:-1,
        emptyITensor(),
        emptyITensor(),
    )
end

function Base.copy(P::ProjNHMPS)
    return ProjNHMPS(
        copy(P.projL), copy(P.projR), P.cached_site_range, copy(P.envL), copy(P.envR)
    )
end

ITensorMPS.nsite(P::ProjNHMPS) = ITensorMPS.nsite(P.projL)

ITensorMPS.site_range(P::ProjNHMPS) = ITensorMPS.site_range(P.projL)

function ITensorMPS.set_nsite!(P::ProjNHMPS, nsite)
    ITensorMPS.set_nsite!(P.projL, nsite)
    ITensorMPS.set_nsite!(P.projR, nsite)
    return P
end

Base.length(P::ProjNHMPS) = length(P.projL)

function getLR(P::ITensorMPS.ProjMPS)
    L = dag(prime(P.M[P.lpos + 1], "Link"))
    !isnothing(lproj(P)) && (L *= lproj(P))

    if ITensorMPS.nsite(P) == 1
        R = ITensor(1)
        !isnothing(rproj(P)) && (R *= rproj(P))
    else
        R = dag(prime(P.M[P.rpos - 1], "Link"))
        !isnothing(rproj(P)) && (R *= rproj(P))
    end

    return L, R
end

function updateenvironments(P::ProjNHMPS)
    ITensorMPS.site_range(P) == P.cached_site_range && return nothing 
        
    Lr, Rr = getLR(P.projR)
    Ll, Rl = getLR(P.projL)

    P.envR = Lr * Rr
    P.envL = Ll * Rl

    P.cached_site_range = ITensorMPS.site_range(P)

    return nothing
end

function adjointproduct(P::ProjNHMPS, v::ITensor)::ITensor
    if ITensorMPS.nsite(P) != 2
        error("Only two-site ProjMPS currently supported")
    end

    updateenvironments(P)

    pv = scalar(P.envR * v)

    Mv = pv * dag(P.envL)

    return noprime(Mv)
end

function product(P::ProjNHMPS, v::ITensor)::ITensor
    # if ITensorMPS.nsite(P) != 2
    #     error("Only two-site ProjMPS currently supported")
    # end

    updateenvironments(P)

    pv = scalar(P.envL * v)

    Mv = pv * dag(P.envR)

    return noprime(Mv)
end

ITensorMPS.position!(P::ProjNHMPS, psi::MPS, pos::Int) = ITensorMPS.position!(P, psi, psi, pos)

function ITensorMPS.position!(P::ProjNHMPS, psil::MPS, psir::MPS, pos::Int)
    ITensorMPS.position!(P.projR, psil, pos)
    ITensorMPS.position!(P.projL, psir, pos)

    return P
end

function ITensorMPS.checkflux(P::ProjNHMPS)
    checkflux(P.projL)
    checkflux(P.projR)
    return nothing
end