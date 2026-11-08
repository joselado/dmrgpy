# Native single-site VUMPS ground state of a translation-invariant chain, for
# Infinite_Many_Body_Chain(..., itensor_version="julia_live").
#
# A dense-array port of pyitensor/vumps.py's grouped (n_uc in {1,2},
# nearest-cell) path, on plain Julia arrays with LinearAlgebra and the
# KrylovKit. The Hamiltonian arrives as the
# grouped automaton W[a,i,o,b] of pyitensor/idmrg.py's _build_periodic_mpo
# (channel 1 = "S", the identity start, channel 2 = "F", the accumulator,
# channels 3.. = one pending channel per open bond); everything else (the
# environments, the two local eigenproblems, the polar updates, the D-ramp
# with restarts, the static observables) runs here.
#
# Conventions, identical to vumps.py so that the arrays mean the same thing
# on both sides: a site tensor is A[l,s,r]; an operator matrix is M[i,o],
# applied to a ket as out[l,o,r] = sum_i M[i,o] A[l,i,r]; an environment
# matrix is X[ket,bra]. Wrapped in a module because every mpsjulialive file
# is evaluated into Main.

module DmrgpyVumps

using LinearAlgebra
using Random
# KrylovKit is declared in juliapkg.json. It is loaded by its UUID, which
# also finds it when an environment holds it only as a dependency of
# ITensorMPS. If it cannot be loaded, the module still loads: the dense
# branches work and only a matrix-free branch raises, when it is reached.
const KK = try
    Base.require(Base.PkgId(Base.UUID("0b1a1467-8014-51b9-945f-bf0ae24f4b77"), "KrylovKit"))
catch
    nothing
end
krylovkit() = KK === nothing ?
    error("DmrgpyVumps: KrylovKit could not be loaded; the matrix-free " *
          "branch (bond dimension above the dense thresholds) needs it") : KK

const CT = ComplexF64

# Below these sizes a linear map is built densely and solved directly;
# above them KrylovKit runs matrix-free. Refs, so that a test can force the
# matrix-free branch on a small problem and compare the two.
const DENSE_TRANSFER_MAX = Ref(24)   # bond dimension D (transfer maps are D^2 x D^2)
const DENSE_HEFF_MAX = Ref(256)      # dimension of the H_AC / H_C problems

function set_dense_thresholds(transfer::Integer, heff::Integer)
    old = (DENSE_TRANSFER_MAX[], DENSE_HEFF_MAX[])
    DENSE_TRANSFER_MAX[] = transfer
    DENSE_HEFF_MAX[] = heff
    return old
end

# --- primitives --------------------------------------------------------

function apply_op(M::AbstractMatrix, A::Array{CT,3})
    D1, d, D2 = size(A)
    dout = size(M, 2)
    out = zeros(CT, D1, dout, D2)
    @inbounds for i in 1:d, o in 1:dout
        m = M[i, o]
        m == 0 && continue
        @views out[:, o, :] .+= m .* A[:, i, :]
    end
    return out
end

# Y[r,r'] = sum A[l,o,r] X[l,l'] conj(B[l',o,r'])   (X[ket,bra] pushed right)
function transfer_left(X::AbstractMatrix, A::Array{CT,3}, B::Array{CT,3})
    Y = zeros(CT, size(A, 3), size(B, 3))
    @views for o in 1:size(A, 2)
        Y .+= transpose(A[:, o, :]) * X * conj(B[:, o, :])
    end
    return Y
end

# X[l,l'] = sum A[l,o,r] Y[r,r'] conj(B[l',o,r'])   (Y[ket,bra] pushed left)
function transfer_right(Y::AbstractMatrix, A::Array{CT,3}, B::Array{CT,3})
    X = zeros(CT, size(A, 1), size(B, 1))
    @views for o in 1:size(A, 2)
        X .+= A[:, o, :] * Y * adjoint(B[:, o, :])
    end
    return X
end

# vec(transfer_left(X,A,A)) = TL * vec(X), and likewise for transfer_right
function transfer_left_matrix(A::Array{CT,3})
    D = size(A, 1)
    T = zeros(CT, D * D, D * D)
    @views for o in 1:size(A, 2)
        T .+= kron(adjoint(A[:, o, :]), transpose(A[:, o, :]))
    end
    return T
end

function transfer_right_matrix(A::Array{CT,3})
    D = size(A, 1)
    T = zeros(CT, D * D, D * D)
    @views for o in 1:size(A, 2)
        T .+= kron(conj(A[:, o, :]), A[:, o, :])
    end
    return T
end

mul_right(A::Array{CT,3}, C::AbstractMatrix) =
    reshape(reshape(A, :, size(A, 3)) * C, size(A, 1), size(A, 2), size(C, 2))
mul_left(C::AbstractMatrix, A::Array{CT,3}) =
    reshape(C * reshape(A, size(A, 1), :), size(C, 1), size(A, 2), size(A, 3))

function polar_unitary(M::AbstractMatrix)
    F = svd(M)
    return F.U * F.Vt
end

pair(X::AbstractMatrix, Y::AbstractMatrix) = sum(X .* Y)   # <X|Y> over (ket,bra)

# --- fixed points ------------------------------------------------------

# Dominant fixed point of a transfer map, bond candidate first: in mixed
# canonical gauge the candidate is an exact algebraic fixed point, so a
# residual test is a yes/no question and no eigensolve runs at
# convergence. Only when the residual says the gauge relation does not hold
# (early iterations) is the guarded eigensolver used.
function fixed_point(cand::Matrix{CT}, apply, matrix_builder, D::Int)
    tc = tr(cand)
    c = abs(tc) > 1e-300 ? cand / tc : nothing
    if c !== nothing
        c = (c + c') / 2
    end
    residual(x) = norm(apply(x) - x) / max(norm(x), eps())
    res_c = c === nothing ? Inf : residual(c)
    res_c <= 1e-6 && return c
    v = try
        dominant_eigvec(apply, matrix_builder, D)
    catch err
        # the guarded eigensolver refused (a degenerate dominant
        # eigenvalue, or no density-matrix eigenvector): the bond
        # candidate falls in, as in pyitensor/vumps.py
        (err isa ErrorException && occursin("degenerate", err.msg) && c !== nothing) || rethrow()
        return c
    end
    if c !== nothing && res_c <= residual(v)
        return c
    end
    return v
end

function dominant_eigvec(apply, matrix_builder, D::Int)
    if D <= max(DENSE_TRANSFER_MAX[], 2)
        T = matrix_builder()
        F = eigen(T)
        idx = sortperm(abs.(F.values), rev=true)
        if length(idx) > 1 && abs(F.values[idx[2]]) > (1 - 1e-9) * abs(F.values[idx[1]])
            error("DmrgpyVumps: degenerate dominant transfer eigenvalue")
        end
        v = reshape(F.vectors[:, idx[1]], D, D)
    else
        vals, vecs, info = krylovkit().eigsolve(x -> vec(apply(reshape(x, D, D))),
                                    vec(Matrix{CT}(I, D, D)) .+ 0.0im, 2, :LM;
                                    tol=1e-13, krylovdim=min(D * D, 40), maxiter=200)
        if info.converged < 1
            # an unconverged dominant vector is refused like a degenerate
            # one, so fixed_point falls back to the bond candidate, as
            # pyitensor/vumps.py does on ArpackError
            error("DmrgpyVumps: degenerate dominant transfer eigenvalue " *
                  "(the Krylov eigensolver did not converge)")
        end
        if length(vals) > 1 && abs(vals[2]) > (1 - 1e-9) * abs(vals[1])
            error("DmrgpyVumps: degenerate dominant transfer eigenvalue")
        end
        v = reshape(vecs[1], D, D)
    end
    v = v / tr(v)
    v = (v + v') / 2
    if !all(isfinite, v) || norm(v) > 1e13
        error("DmrgpyVumps: degenerate dominant transfer eigenvalue (the fixed " *
              "point is not a density matrix)")
    end
    return v
end

# r: right fixed point of the AL transfer; l: left fixed point of the AR one
function transfer_fixed_points(AL, AR, C)
    D = size(AL, 1)
    r = fixed_point(C * C', Y -> transfer_right(Y, AL, AL),
                    () -> transfer_right_matrix(AL), D)
    l = fixed_point(Matrix(transpose(C' * C)), X -> transfer_left(X, AR, AR),
                    () -> transfer_left_matrix(AR), D)
    return r, l
end

# --- environments ------------------------------------------------------

# W blocks (a, b, M) with a nonzero operator matrix
function w_blocks(W::Array{CT,4})
    Dw = size(W, 1)
    blocks = Tuple{Int,Int,Matrix{CT}}[]
    for a in 1:Dw, b in 1:Dw
        M = W[a, :, :, b]
        if maximum(abs.(M)) > 0
            push!(blocks, (a, b, Matrix(M)))
        end
    end
    return blocks
end

# GMRES Krylov dimension and restarts; a Ref so that a test can starve it
# and check the dense fallback. LINSOLVE_FALLBACKS counts the fallbacks.
const LINSOLVE_LIMITS = Ref((60, 200))
const LINSOLVE_FALLBACKS = Ref(0)

# solve X - T(X) + I*<X|fp> = rhs
function regularized_solve(rhs::Matrix{CT}, apply, matrix_builder, fp, D::Int)
    Id = Matrix{CT}(I, D, D)
    if D <= max(DENSE_TRANSFER_MAX[], 2)
        T = matrix_builder()
        A = Matrix{CT}(I, D * D, D * D) - T + vec(Id) * transpose(vec(fp))
        return reshape(A \ vec(rhs), D, D)
    end
    f = x -> begin
        X = reshape(x, D, D)
        vec(X - apply(X) + Id * pair(X, fp))
    end
    kdim, kiter = LINSOLVE_LIMITS[]
    sol, info = krylovkit().linsolve(f, vec(rhs), vec(rhs); tol=1e-14 * max(norm(rhs), 1e-300),
                            krylovdim=min(D * D, kdim), maxiter=kiter)
    LINSOLVE_FALLBACKS[] += info.converged < 1
    if info.converged < 1
        # GMRES did not converge: solve the same system densely, as
        # idmrg_excitations._solve_linear_map does, rather than use a
        # wrong environment
        T = matrix_builder()
        A = Matrix{CT}(I, D * D, D * D) - T + vec(Id) * transpose(vec(fp))
        return reshape(A \ vec(rhs), D, D)
    end
    return reshape(sol, D, D)
end

# Left environments L[a] (bond after a left half of AL's, in channel a) and
# right environments R[a] (bond before a right half of AR's, in channel a);
# reach-1 automaton: pending channels are entered only from S and left only
# into F. e_L / e_R: the energy per (grouped) site seen from each side.
function environments(AL, AR, C, blocks, Dw::Int)
    D = size(AL, 1)
    r, l = transfer_fixed_points(AL, AR, C)
    Id = Matrix{CT}(I, D, D)
    L = [zeros(CT, D, D) for _ in 1:Dw]
    R = [zeros(CT, D, D) for _ in 1:Dw]
    L[1] = copy(Id)
    R[2] = copy(Id)
    for (a, b, M) in blocks
        if a == 1 && b >= 3
            L[b] .+= transfer_left(L[1], apply_op(M, AL), AL)
        end
        if b == 2 && a >= 3
            R[a] .+= transfer_right(R[2], apply_op(M, AR), AR)
        end
    end
    YL = zeros(CT, D, D)
    YR = zeros(CT, D, D)
    for (a, b, M) in blocks
        if b == 2 && a != 2
            YL .+= transfer_left(L[a], apply_op(M, AL), AL)
        end
        if a == 1 && b != 1
            YR .+= transfer_right(R[b], apply_op(M, AR), AR)
        end
    end
    eL = real(pair(YL, r))
    eR = real(pair(l, YR))
    L[2] = regularized_solve(YL - eL * Id, X -> transfer_left(X, AL, AL),
                             () -> transfer_left_matrix(AL), r, D)
    R[1] = regularized_solve(YR - eR * Id, Y -> transfer_right(Y, AR, AR),
                             () -> transfer_right_matrix(AR), l, D)
    L[2] = (L[2] + L[2]') / 2
    R[1] = (R[1] + R[1]') / 2
    return L, R, eL, eR
end

# --- effective Hamiltonians -------------------------------------------

function h_ac(X::Array{CT,3}, L, R, blocks)
    D1, d, D2 = size(X)
    out = zeros(CT, D1, d, D2)
    for (a, b, M) in blocks
        Xo = apply_op(M, X)
        T1 = transpose(L[a]) * reshape(Xo, D1, :)
        out .+= reshape(reshape(T1, :, D2) * R[b], D1, d, D2)
    end
    return out
end

function h_c(X::AbstractMatrix, L, R)
    out = zeros(CT, size(X))
    for a in eachindex(L)
        out .+= transpose(L[a]) * X * R[a]
    end
    return out
end

# Lowest eigenvector of a Hermitian linear map on arrays of shape `shp`,
# warm-started from x0, stopped on the Ritz residual. Dense below
# DENSE_HEFF_MAX.
function lowest_eigvec(apply, x0::Array{CT}, rtol::Float64, niter::Int)
    shp = size(x0)
    n = length(x0)
    if n <= max(DENSE_HEFF_MAX[], 4)
        H = zeros(CT, n, n)
        e = zeros(CT, n)
        for j in 1:n
            e .= 0
            e[j] = 1
            H[:, j] = vec(apply(reshape(e, shp)))
        end
        F = eigen(Hermitian((H + H') / 2))
        return F.values[1], reshape(F.vectors[:, 1], shp)
    end
    v0 = vec(x0) / norm(x0)
    vals, vecs, _ = krylovkit().eigsolve(x -> vec(apply(reshape(x, shp))), v0, 1, :SR;
                                ishermitian=true, tol=rtol,
                                krylovdim=min(n, max(niter, 20)), maxiter=100)
    return real(vals[1]), reshape(vecs[1] / norm(vecs[1]), shp)
end

# remove the whole relative phase against the previous vector
function align_phase(new::Array{CT}, old::Array{CT})
    ov = dot(old, new)
    return abs(ov) > 0 ? new * conj(ov / abs(ov)) : new
end

function update_AL_AR(AC::Array{CT,3}, C::Matrix{CT})
    D, d, _ = size(AC)
    Ul = polar_unitary(reshape(AC, D * d, D))
    Uc = polar_unitary(C)
    AL = reshape(Ul * Uc', D, d, D)
    Ur = polar_unitary(reshape(AC, D, d * D))
    AR = reshape(Uc' * Ur, D, d, D)
    return AL, AR
end

function gauge_mismatch(AC, C, AL, AR)
    n = norm(AC)
    n == 0 && return 0.0
    return (norm(AC - mul_right(AL, C)) + norm(AC - mul_left(C, AR))) / n
end

# --- initial states ------------------------------------------------------

function random_state(D::Int, d::Int, rng)
    AL = reshape(polar_unitary(randn(rng, CT, D * d, D)), D, d, D)
    AR = reshape(polar_unitary(randn(rng, CT, D, d * D)), D, d, D)
    C = Matrix{CT}(I, D, D) / sqrt(D)
    return AL, AR, C
end

# the previous rung's state embedded in the top-left block plus noise
function grown_state(D::Int, d::Int, prev, rng; noise=0.05)
    Do = size(prev.AL, 1)
    AL0 = noise * randn(rng, CT, D, d, D)
    AR0 = noise * randn(rng, CT, D, d, D)
    AL0[1:Do, :, 1:Do] .+= prev.AL
    AR0[1:Do, :, 1:Do] .+= prev.AR
    AL = reshape(polar_unitary(reshape(AL0, D * d, D)), D, d, D)
    AR = reshape(polar_unitary(reshape(AR0, D, d * D)), D, d, D)
    C = noise * randn(rng, CT, D, D) / D
    C[1:Do, 1:Do] .+= prev.C
    C ./= norm(C)
    return AL, AR, C
end

# --- the solver ----------------------------------------------------------

mutable struct VumpsState
    AL::Array{CT,3}
    AR::Array{CT,3}
    C::Matrix{CT}
    AC::Array{CT,3}
    e0::Float64          # energy per physical site
    converged::Bool
    mismatch::Float64
    niter::Int
    D::Int
    d::Int
end

function single_run(blocks, Dw, n_uc, D, d, tol, maxiter, niter, hunit, init, verbose)
    AL, AR, C = init
    AC = mul_right(AL, C)
    AC ./= norm(AC)
    C ./= norm(C)
    converged = false
    mismatch = Inf
    lam = 1.0
    it = 0
    for outer it in 1:maxiter
        L, R, eL, eR = environments(AL, AR, C, blocks, Dw)
        rtol = (tol / 10) * min(hunit, max(1.0, abs(lam)))
        lam, ACn = lowest_eigvec(X -> h_ac(X, L, R, blocks), AC, rtol, niter)
        _, Cn = lowest_eigvec(X -> h_c(X, L, R), C, rtol, niter)
        ACn = align_phase(ACn, AC)
        Cn = align_phase(Cn, C)
        mismatch = gauge_mismatch(ACn, Cn, AL, AR)
        if verbose
            println("julia vumps D=$D iter $it: e0=$((eL + eR) / 2 / n_uc) mismatch=$mismatch")
        end
        AC, C = ACn, Cn
        AL, AR = update_AL_AR(AC, C)
        if mismatch < tol
            converged = true
            break
        end
    end
    _, _, eL, eR = environments(AL, AR, C, blocks, Dw)
    return VumpsState(AL, AR, C, AC, (eL + eR) / 2 / n_uc, converged, mismatch, it, D, d)
end

function d_ramp(D::Int)
    out = Int[]
    x = 1
    while x < D
        push!(out, x)
        x *= 2
    end
    push!(out, D)
    return out
end

better(a, b) = a === nothing ? b : b === nothing ? a :
    (a.converged != b.converged ? (a.converged ? a : b) : (a.e0 < b.e0 ? a : b))

"""
    vumps_ground_state(W, n_uc, D, tol, maxiter, niter, nrestarts, hunit, verbose)

W: grouped automaton W[a,i,o,b] (channel 1 = S, 2 = F, 3.. pending,
reach 1). Ramps D = 1, 2, 4, ..., D, each rung warm-started from the
previous best (noise padding), `nrestarts` attempts at the target D (at most
3 on the way), keeping converged over unconverged and then the lowest energy,
with the variational safety net of pyitensor/vumps.py: a rung whose best
energy is above an already-known smaller-D energy gets up to 2*nrestarts
more attempts from the same anchor.
"""
function vumps_ground_state(W_in, n_uc::Integer, D::Integer, tol::Real,
                            maxiter::Integer, niter::Integer, nrestarts::Integer,
                            hunit::Real, verbose::Bool)
    W = Array{CT,4}(W_in)
    Dw = size(W, 1)
    d = size(W, 2)
    blocks = w_blocks(W)
    for (a, b, _) in blocks
        if a >= 3 && b >= 3
            error("DmrgpyVumps: automaton is not reach-1 (pending-to-pending block)")
        end
    end
    rng = Random.default_rng()
    function attempt(Dc, init)
        try
            return single_run(blocks, Dw, Int(n_uc), Dc, d, Float64(tol), Int(maxiter),
                              Int(niter), Float64(hunit), init, verbose)
        catch err
            # only the guarded fixed-point eigensolver's refusal skips an
            # attempt; anything else is a bug and propagates
            if err isa ErrorException && occursin("degenerate", err.msg)
                verbose && println("julia vumps D=$Dc attempt failed: $err")
                return nothing
            end
            rethrow()
        end
    end
    best = nothing
    prev = nothing
    best_e = nothing
    margin = 1e-6 * min(1.0, Float64(hunit))
    for Dc in d_ramp(Int(D))
        nhere = Dc == D ? nrestarts : min(nrestarts, 3)
        local_best = nothing
        for k in 1:nhere
            init = (k == 1 && prev !== nothing) ? grown_state(Dc, d, prev, rng) :
                   random_state(Dc, d, rng)
            local_best = better(attempt(Dc, init), local_best)
        end
        if local_best === nothing
            error("DmrgpyVumps: every attempt at D=$Dc failed (degenerate " *
                  "transfer-matrix spectrum); try increasing nrestarts")
        end
        if best_e !== nothing && local_best.e0 > best_e + margin && prev !== nothing
            for _ in 1:(2 * nrestarts)
                local_best.e0 <= best_e + margin && break
                local_best = better(attempt(Dc, grown_state(Dc, d, prev, rng)), local_best)
            end
        end
        prev = local_best
        best_e = best_e === nothing ? local_best.e0 : min(best_e, local_best.e0)
        best = local_best
    end
    return best
end

# --- static observables ----------------------------------------------------

# <M> on the grouped site, from AC (exact single-site reduced state)
function onsite_expectation(st::VumpsState, M_in)
    M = Matrix{CT}(M_in)
    AC = st.AC
    return dot(AC, apply_op(M, AC)) / real(dot(AC, AC))
end

# <Mi (cell 0) ... Mj (cell `offset`)>, offset >= 1: AC carries the first
# cell, AR every cell to its right (exactly right-orthonormal, so the open
# right bond closes with a plain trace). Mcross is the Jordan-Wigner string
# on the cells strictly between, `nothing` for bosonic operators.
function two_point_correlator(st::VumpsState, Mi_in, Mcross_in, Mj_in, offset::Integer)
    AC, AR = st.AC, st.AR
    D = size(AC, 1)
    Id = Matrix{CT}(I, D, D)
    X = transfer_left(Id, apply_op(Matrix{CT}(Mi_in), AC), AC) / real(dot(AC, AC))
    ARx = Mcross_in === nothing ? AR : apply_op(Matrix{CT}(Mcross_in), AR)
    for _ in 1:(offset - 1)
        X = transfer_left(X, ARx, AR)
    end
    X = transfer_left(X, apply_op(Matrix{CT}(Mj_in), AR), AR)
    return tr(X)
end

# diagnostics: residuals of the two fixed-point candidates
function candidate_residuals(st::VumpsState)
    AL, AR, C = st.AL, st.AR, st.C
    r = C * C'
    r /= tr(r)
    l = Matrix(transpose(C' * C))
    l /= tr(l)
    lw = C' * C
    lw /= tr(lw)
    return (norm(transfer_right(r, AL, AL) - r), norm(transfer_left(l, AR, AR) - l),
            norm(transfer_left(lw, AR, AR) - lw))
end

end # module DmrgpyVumps
