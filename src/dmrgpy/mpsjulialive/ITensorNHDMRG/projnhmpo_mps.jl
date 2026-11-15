using ITensors

mutable struct ProjNHMPO_MPS
    PH::ProjNHMPO
    pm::Vector{ProjNHMPS}
    weight::Float64
end

Base.copy(P::ProjNHMPO_MPS) = ProjNHMPO_MPS(copy(P.PH), copy.(P.pm), P.weight)

function ProjNHMPO_MPS(H::MPO, mpsvl::Vector{MPS}, mpsvr::Vector{MPS}; weight=1.0, keepadj=true)
    @assert length(mpsvl) == length(mpsvr)
    return ProjNHMPO_MPS(
        ProjNHMPO(H; keepadj), [ProjNHMPS(mpsvl[i], mpsvr[i]) for i in eachindex(mpsvl)], weight
    )
end

(P::ProjNHMPO_MPS)(v::ITensor) = product(P, v)

ITensorMPS.nsite(P::ProjNHMPO_MPS) = ITensorMPS.nsite(P.PH)

function ITensorMPS.set_nsite!(Ps::ProjNHMPO_MPS, nsite)
    ITensorMPS.set_nsite!(Ps.PH, nsite)
    for P in Ps.pm
        ITensorMPS.set_nsite!(P, nsite)
    end
    return Ps
end

Base.length(P::ProjNHMPO_MPS) = length(P.PH)

function ITensorMPS.site_range(P::ProjNHMPO_MPS)
    r = ITensorMPS.site_range(P.PH)
    @assert all(m -> ITensorMPS.site_range(m) == r, P.pm)
    return r
end

function adjointproduct(P::ProjNHMPO_MPS, v::ITensor)::ITensor
    Pv = adjointproduct(P.PH, v)
    for p in P.pm
        Pv += P.weight * adjointproduct(p, v)
    end
    return Pv
end

function product(P::ProjNHMPO_MPS, v::ITensor)::ITensor
    Pv = product(P.PH, v)
    for p in P.pm
        Pv += P.weight * product(p, v)
    end
    return Pv
end

function Base.eltype(P::ProjNHMPO_MPS)
    elT = eltype(P.PH)
    for p in P.pm
        elT = promote_type(elT, eltype(p))
    end
    return elT
end

Base.size(P::ProjNHMPO_MPS) = size(P.H)

function ITensorMPS.position!(P::ProjNHMPO_MPS, psil::MPS,psir::MPS, pos::Int)
    ITensorMPS.position!(P.PH, psil, psir, pos)
    for p in P.pm
        ITensorMPS.position!(p, psil, psir, pos)
    end
    return P
end

ITensorMPS.position!(P::ProjNHMPO_MPS, psi::MPS, pos::Int) = ITensorMPS.position!(P, psi, psi, pos)

ITensorMPS.noiseterm(P::ProjNHMPO_MPS, thetal::ITensor, thetar::ITensor, ortho::String) = ITensorMPS.noiseterm(P.PH, thetal, thetar, ortho)
ITensorMPS.noiseterm(P::ProjNHMPO_MPS, theta::ITensor, ortho::String) = ITensorMPS.noiseterm(P.PH, theta, theta, ortho)

function ITensorMPS.checkflux(P::ProjNHMPO_MPS)
    checkflux(P.PH)
    foreach(checkflux, P.pm)
    return nothing
end