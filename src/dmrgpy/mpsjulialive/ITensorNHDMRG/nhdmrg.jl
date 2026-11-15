using ITensors: @timeit_debug, @debug_check, scalartype, Printf, Algorithm, @Algorithm_str
using ITensorMPS: check_hascommoninds, ProjMPO, position!
using Printf: @printf
using LinearAlgebra

nhdmrg(H::MPO, psi::MPS, sweeps::Sweeps; kwargs...) = nhdmrg(H, psi, psi, sweeps; kwargs...)

function nhdmrg(H::MPO, psil0::MPS, psir0::MPS, sweeps::Sweeps; kwargs...)
    check_hascommoninds(siteinds, H, psir0)
    check_hascommoninds(siteinds, H, psil0')
    check_hascommoninds(siteinds, psir0, psil0)
    # Permute the indices to have a better memory layout
    # and minimize permutations
    H = permute(H, (linkind, siteinds, linkind))
    PH = ProjNHMPO(H)
    
    return nhdmrg(PH, psil0, psir0, sweeps; kwargs...)
end

function nhdmrg(
    H::MPO,
    Msl::Vector{MPS},
    Msr::Vector{MPS},
    psi::MPS,
    sweeps::Sweeps;
    weight=true,
    kwargs...,
)
    return nhdmrg(H, Msl, Msr, psi, psi, sweeps; weight, kwargs...)
end

function nhdmrg(
    H::MPO,
    Msl::Vector{MPS},
    Msr::Vector{MPS},
    psil0::MPS,
    psir0::MPS,
    sweeps::Sweeps;
    weight=true,
    kwargs...,
)
    check_hascommoninds(siteinds, H, psir0)
    check_hascommoninds(siteinds, H, psil0')
    check_hascommoninds(siteinds, psir0, psil0)
    # Permute the indices to have a better memory layout
    # and minimize permutations
    H = permute(H, (linkind, siteinds, linkind))
    PH = ProjNHMPO_MPS(H, Msl, Msr; weight)
    return nhdmrg(PH, psil0, psir0, sweeps; kwargs...)
end

"""
    nhdmrg(H::MPO, psir0::MPS, psil0::MPS, sweeps::Sweeps; kwargs...)

Use a non-Hermitian variant of the density matrix renormalization group 
(NH-DMRG) algorithm to optimize the matrix product states (MPS) 
such that they are the left- and right-eigenvector of lowest eigenvalue 
of a non-Hermitian matrix `H`, represented as a matrix product operator (MPO).

    nhdmrg(H::MPO, Msl::Vector{MPS}, Msr::Vector{MPS}, psir0::MPS, psil0::MPS, sweeps::Sweeps; weight=1.0, kwargs...)
    
Use a non-Hermitian variant of the density matrix renormalization group 
(NH-DMRG) algorithm to optimize the matrix product states (MPS) 
such that they are the left- and right-eigenvector of lowest eigenvalue 
of a non-Hermitian matrix `H`, subject to the constraint that the MPS
are orthogonal to each of the left- and right-MPS provided in the Vector 
`Msl` and `Msr`. The orthogonality constraint is approximately enforced by
adding to `H` terms of the form `w|Mr1><Ml1| + w|Mr2><Ml2| + ...` where 
`Mls=[Ml1, Ml2, ...]`, `Mrs=[Mr1, Mr2, ...]` and `w` is the "weight" parameter, 
which can be adjusted through the optional `weight` keyword argument.

!!! note
    `nhdmrg` will report the energy of the operator
    `H + w|Mr1><Ml1| + w|Mr2><Ml2| + ...`, not the operator `H`.
    If you want the expectation value of the MPS eigenstate
    with respect to just `H`, you can compute it yourself with
    an observer or after DMRG is run with `inner(psil', H, psir)`.

The MPS `psil0` and `psir0` is used to initialize the MPS to be optimized.

The number of sweeps and accuracy parameters can be passed through a 
`Sweeps` object.

The NH-DMRG algorithm allows choosing both the eigenvalue solver by passing 
the `alg` keyword and the truncation algorithm by passing the `biorthoalg`
keyword. The eigenvalue algorithm `alg` may be `onesided`, such that the for each 
block the eigenvalue problem ``A |x> = λ |x>`` and ``A† |y> = λ* |y>`` are 
solved, or `twosided`, such that ``<y| A |x> = λ <y|x>`` are solved simultaneously
using a two-sided Krylov approach [1]. The truncation algorithm `biorthoalg` currently 
supports the biorthogonal block method `biorthoblock` [2] and the left-right density 
matrix method `lrdensity` [3].

[1] https://doi.org/10.1137/16M1078987
[2] https://doi.org/10.48550/arXiv.2401.15000
[3] https://doi.org/10.1103/PhysRevB.105.205125

Returns:

  - `energy::Number` - eigenvalue of the optimized MPS
  - `psil::MPS` - optimized left MPS
  - `psir::MPS` - optimized right MPS

Optional keyword arguments:

  - `alg::String = "onesided"` - local eigenvalue algorithm, either `"onesided"` or `"twosided"`
  - `biorthoalg::String = "fidelity"` - orthogonalization algorithm, eiter `"biorthoblock"` 
     or `"fidelity"`
  - `isbiortho::Bool=false` - if `true` the input MPS are not biorthogonalized before starting 
     the sweeps 
  - `eigsolve_krylovdim::Int = 3` - maximum dimension of Krylov space used to
     locally solve the eigenvalue problem. Try setting to a higher value if
     convergence is slow or the Hamiltonian is close to a critical point. [^krylovkit]
  - `eigsolve_tol::Number = 1e-14` - Krylov eigensolver tolerance. [^krylovkit]
  - `eigsolve_maxiter::Int = 1` - number of times the Krylov subspace can be
     rebuilt. [^krylovkit]
  - `eigsolve_verbosity::Int = 0` - verbosity level of the Krylov solver.
     Warning: enabling this will lead to a lot of outputs to the terminal. [^krylovkit]
  - `outputlevel::Int = 1` - larger outputlevel values make DMRG print more
     information and 0 means no output.
  - `observer` - object implementing the [Observer](@ref observer) interface
     which can perform measurements and stop DMRG early.
  - `biorthokwargs` - additional kwargs for the biorthogonalization routine, 
     see `decomposerho`.
"""
function nhdmrg(
    PH::Union{ProjNHMPO,ProjNHMPO_MPS,ProjNHMPS},
    psil0::MPS,
    psir0::MPS,
    sweeps::Sweeps;
    alg="onesided",
    biorthoalg="fidelity",
    isbiortho::Bool=false,
    observer=NoObserver(),
    outputlevel=1,
    # eigsolve kwargs
    eigsolve_tol=1e-14,
    eigsolve_krylovdim=3,
    eigsolve_maxiter=1,
    eigsolve_verbosity=0,
    eigsolve_which_eigenvalue=:SR,
    biorthokwargs...,
)
    if length(psir0) == 1
        error(
            "`dmrg` currently does not support system sizes of 1. You can diagonalize the MPO tensor directly with tools like `LinearAlgebra.eigen`, `KrylovKit.eigsolve`, etc.",
        )
    end

    @debug_check begin
        # Debug level checks
        # Enable with ITensors.enable_debug_checks()
        checkflux(psir0)
        checkflux(PH)
    end

    psir = deepcopy(psir0)
    psil = deepcopy(psil0)
    N = length(psir0)

    if !isbiortho
        psil, psir = biorthogonalize!(
            psil,
            psir,
            biorthoalg;
            maxdim=maxdim(sweeps, 1),
            mindim=mindim(sweeps, 1),
            cutoff=cutoff(sweeps, 1),
            biorthokwargs...,
        )
    end

    overlap = inner(psil, psir)
    if isreal(overlap)
        psil /= sqrt(abs(overlap))
        psir /= sqrt(abs(overlap)) * sign(overlap)   
    else
        psil /= sqrt(conj(overlap))
        psir /= sqrt(overlap)
    end

    PH = ITensorMPS.position!(PH, psil, psir, 1)
    
    energy = 0.0

    for sw in 1:nsweep(sweeps)
        sw_time = @elapsed begin
            maxtruncerr = 0.0

            for (b, ha) in sweepnext(N)
                @debug_check begin
                    checkflux(psir)
                    checkflux(PH)
                end

                PH = ITensorMPS.position!(PH, psil, psir, b)

                @debug_check begin
                    checkflux(psir)
                    checkflux(PH)
                end

                begin 
                    Θl = psil[b] * psil[b + 1]
                    Θr = psir[b] * psir[b + 1]
                end 

                energy, nΘl, nΘr = nhproblemsolver!(
                    Algorithm(alg),
                    PH,
                    Θl,
                    Θr;
                    eigsolve_tol,
                    eigsolve_krylovdim,
                    eigsolve_maxiter,
                    eigsolve_verbosity,
                    eigsolve_which_eigenvalue
                )

                nΘl /= sqrt(inner(nΘl, nΘl))
                nΘr /= sqrt(inner(nΘr, nΘr))

                begin 
                    overlap = inner(nΘl, nΘr)
                    if isreal(overlap)
                        nΘl /= sqrt(abs(overlap))
                        nΘr /= sqrt(abs(overlap)) * sign(overlap)
                    else
                        nΘl /= sqrt(conj(overlap))
                        nΘr /= sqrt(overlap)
                    end
                end

                ortho = ha == 1 ? "left" : "right"

                drho = nothing
                if noise(sweeps, sw) > 0
                    # Use noise term when determining new MPS basis.
                    # This is used to preserve the element type of the MPS.
                    elt = real(scalartype(psir))
                    drho = elt(noise(sweeps, sw)) * noiseterm(PH, nΘl, nΘr, ortho)
                end

                spec = nhreplacebond!(
                    psil,
                    psir,
                    b,
                    nΘl, 
                    nΘr,
                    biorthoalg;
                    ortho,
                    eigen_perturbation=drho,
                    maxdim=maxdim(sweeps, sw),
                    mindim=mindim(sweeps, sw),
                    cutoff=cutoff(sweeps, sw),
                    biorthokwargs...,
                )

                maxtruncerr = max(maxtruncerr, spec.truncerr)

                @debug_check begin
                    checkflux(psir)
                    checkflux(PH)
                end

                if outputlevel >= 2
                    @printf(
                        "Sweep %d, half %d, bond (%d,%d) energy=%s\n",
                        sw,
                        ha,
                        b,
                        b + 1,
                        energy
                    )
                    @printf(
                        "  Truncated using cutoff=%.1E maxdim=%d mindim=%d\n",
                        cutoff(sweeps, sw),
                        maxdim(sweeps, sw),
                        mindim(sweeps, sw)
                    )
                    flush(stdout)
                end

                sweep_is_done = (b == 1 && ha == 2)
                measure!(
                    observer;
                    energy,
                    psi=(psil, psir),
                    projected_operator=PH,
                    bond=b,
                    sweep=sw,
                    half_sweep=ha,
                    spec,
                    outputlevel,
                    sweep_is_done,
                )
            end
        end
        if outputlevel >= 1
            @printf(
                "After sweep %d energy=%s  maxlinkdim=%d maxerr=%.2E time=%.3f\n",
                sw,
                energy,
                maxlinkdim(psil),
                maxtruncerr,
                sw_time
            )
            flush(stdout)
        end
        isdone = checkdone!(observer; energy, psir, sweep=sw, outputlevel)
        isdone && break
    end
    return (energy, psil, psir)
end

function nhdmrg(
    x1,
    x2,
    x3,
    psil0::MPS,
    psir0::MPS;
    nsweeps,
    maxdim=ITensorMPS.default_maxdim(),
    mindim=ITensorMPS.default_mindim(),
    cutoff=ITensorMPS.default_cutoff(Float64),
    noise=ITensorMPS.default_noise(),
    kwargs...,
)
    return nhdmrg(
        x1,
        x2,
        x3,
        psil0,
        psir0,
        ITensorMPS._dmrg_sweeps(; nsweeps, maxdim, mindim, cutoff, noise);
        kwargs...,
    )
end

function nhdmrg(
    x1,
    x2,
    x3,
    psi::MPS,
    nsweeps,
    maxdim=ITensorMPS.default_maxdim(),
    mindim=ITensorMPS.default_mindim(),
    cutoff=ITensorMPS.default_cutoff(Float64),
    noise=ITensorMPS.default_noise(),
    kwargs...,
)
    return nhdmrg(
        x1,
        x2,
        x3,
        psi,
        ITensorMPS._dmrg_sweeps(; nsweeps, maxdim, mindim, cutoff, noise);
        kwargs...,
    )
end

function nhdmrg(
    x1,
    psil0::MPS,
    psir0::MPS;
    nsweeps,
    maxdim=ITensorMPS.default_maxdim(),
    mindim=ITensorMPS.default_mindim(),
    cutoff=ITensorMPS.default_cutoff(Float64),
    noise=ITensorMPS.default_noise(),
    kwargs...,
)
    return nhdmrg(
        x1,
        psil0,
        psir0,
        ITensorMPS._dmrg_sweeps(; nsweeps, maxdim, mindim, cutoff, noise);
        kwargs...,
    )
end


function nhdmrg(
    x1,
    psi::MPS;
    nsweeps,
    maxdim=ITensorMPS.default_maxdim(),
    mindim=ITensorMPS.default_mindim(),
    cutoff=ITensorMPS.default_cutoff(Float64),
    noise=ITensorMPS.default_noise(),
    kwargs...,
)
    return nhdmrg(
        x1,
        psi,
        ITensorMPS._dmrg_sweeps(; nsweeps, maxdim, mindim, cutoff, noise);
        kwargs...,
    )
end