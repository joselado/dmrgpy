# Controlled bond expansion (CBE) for one-site DMRG on julia_live.
#
# A port of pyitensor/cbe.py's dmrg_cbe (shrewd selection only), which has
# the references (Gleis, Li, von Delft, PRL 130, 246402 (2023),
# arXiv:2207.14712) and the derivation; the kernel below follows
# expand_left_isometry line by line. One-site DMRG never leaves the space
# set by the current bonds. Before each one-site update the bond to the
# next site is expanded by the few directions of the two-site discarded
# space on which H|psi> has weight, with zero weight in the state, so
# the state is unchanged and the one-site update can use them; the SVD
# that moves the center trims the bond back to maxdim and the cutoff.
#
# Environments are plain ITensor contractions,
#   L_i = L_{i-1} * psi[i] * H[i] * dag(prime(psi[i])),
# kept by hand rather than through ProjMPO, since the expansion rewrites
# psi[i+1] behind any cache. The expansion itself is dense linear algebra
# on Julia arrays, so it needs dense (non-QN) indices, which is what every
# julia_live chain has.
#
# KrylovKit is not a direct dependency of the julia_live project, only of
# ITensorMPS, so it is reached through that module.

import LinearAlgebra
import Random

const _CBE_KK = ITensorMPS.KrylovKit

# Singular values below this fraction of the scale of the object they come
# from are numerical zeros (pyitensor/cbe.py's _NUMERICAL_ZERO).
const _CBE_NUMERICAL_ZERO = 1e-12

_cbe_dim(x) = x === nothing ? 1 : dim(x)
_cbe_prime(x) = x === nothing ? nothing : prime(x)

# T as a dense array in the order of inds (nothing entries skipped), then
# reshaped to dims, which counts a skipped index as 1.
function _cbe_array(T, inds, dims)
  present = [ix for ix in inds if ix !== nothing]
  arr = array(T, present...)
  return reshape(arr, dims...)
end

function _cbe_itensor(arr, inds)
  present = [ix for ix in inds if ix !== nothing]
  return ITensor(reshape(arr, Tuple(dim(ix) for ix in present)), present...)
end

# How many of the descending singular values s to keep: those above
# max(eps, numerical zero)*scale, at most cap (nothing for no cap).
function _cbe_count(s, scale, eps, cap)
  (scale <= 0 || isempty(s)) && return 0
  floor = max(eps, _CBE_NUMERICAL_ZERO)*scale
  keep = count(x -> real(x) > floor, s)
  cap === nothing || (keep = min(keep, cap))
  return max(keep, 0)
end

# Orthonormal rows spanning the row space of C (k, n), as an (r, n) matrix.
function _cbe_row_basis(C)
  F = LinearAlgebra.svd(C)
  r = _cbe_count(F.S, isempty(F.S) ? 0.0 : F.S[1], 0.0, nothing)
  return F.Vt[1:r, :]
end

# U (n, m) with the span of A's columns removed, made orthonormal again,
# numerically zero columns dropped; applied twice (the paper's step 14).
function _cbe_complement(U, A)
  for _ in 1:2
    U = U - A*(A'*U)
  end
  size(U, 2) == 0 && return U
  F = LinearAlgebra.svd(U)
  r = _cbe_count(F.S, 1.0, 0.0, nothing)
  return F.U[:, 1:r]
end

# The shrewd selection: new columns At (nl, D~) for the isometry Am (nl, k)
# of a bond whose center Cm (k, nr) sits on the other side. X3 (nl, dv, k)
# is the half of H|psi> on Am's side, open at the central MPO bond dv and
# at the ket bond k; Y3 (k, dv, nr) is the half on the center's side.
# expand_left_isometry in pyitensor/cbe.py, with numpy's einsum and
# reshapes written as matrix products of the same objects.
function _cbe_kernel(X3, Y3, Am, Cm, dtilde, dpre, eps_pre, eps_final)
  nl, dv, k = size(X3)
  nr = size(Y3, 3)
  empty = zeros(promote_type(eltype(X3), eltype(Y3)), nl, 0)
  dtilde <= 0 && return empty
  Vb = _cbe_row_basis(Cm)                          # (r, nr)
  Xm = reshape(X3, nl, dv*k)                       # columns (v, k)
  Yq = reshape(Y3, k*dv, nr)                       # rows (k, v)
  # (a) the center half projected on (1 - B^dag B), SVD on the ket bond
  Yo = Yq - (Yq*Vb')*Vb
  F1 = LinearAlgebra.svd(reshape(Yo, k, dv*nr))
  r1 = _cbe_count(F1.S, LinearAlgebra.norm(Y3), 0.0, nothing)
  r1 == 0 && return empty
  US = F1.U[:, 1:r1] .* transpose(F1.S[1:r1])     # (k, r1)
  # the other half with these weights on the central MPS bond, projected
  # on (1 - A A^dag)
  Lt = reshape(reshape(X3, nl*dv, k)*US, nl, dv*r1)
  scale_pre = LinearAlgebra.norm(Lt)
  Lo = Lt - Am*(Am'*Lt)
  # (b) truncate the central MPS bond r1 -> D'
  F2 = LinearAlgebra.svd(reshape(Lo, nl*dv, r1))
  r2 = _cbe_count(F2.S, scale_pre, eps_pre, dpre)
  r2 == 0 && return empty
  us = reshape(F2.U[:, 1:r2] .* transpose(F2.S[1:r2]), nl, dv*r2)
  # (c) regroup the MPO bond with the D' kept states: the preselected
  # complement A^pr
  F3 = LinearAlgebra.svd(us)
  r3 = _cbe_count(F3.S, isempty(F3.S) ? 0.0 : F3.S[1], 0.0, nothing)
  Apr = _cbe_complement(F3.U[:, 1:r3], Am)
  size(Apr, 2) == 0 && return empty
  # (d) final selection: H|psi> on A^pr x (1 - B^dag B), MPO bond closed
  Lpr = Apr'*Xm                                    # (x, (v, k))
  nx = size(Lpr, 1)
  Ct = reshape(permutedims(reshape(Lpr, nx, dv, k), (1, 3, 2)), nx, k*dv)*Yq
  scale_fin = LinearAlgebra.norm(Ct)
  Co = Ct - (Ct*Vb')*Vb
  F4 = LinearAlgebra.svd(Co)
  r4 = _cbe_count(F4.S, scale_fin, eps_final, dtilde)
  return Apr*F4.U[:, 1:r4]
end

# The neighbour link of site j on the side away from site o, or nothing at
# a chain end.
function _cbe_outer_link(psi, j, o)
  q = j + (j - o)
  (q < 1 || q > length(psi)) && return nothing
  return commonind(psi[j], psi[q])
end

# Expand the bond between site p (an isometry, the side that grows) and
# site c (the orthogonality center), its neighbour. Ep is the environment
# beyond p, Ec the one beyond c. Rewrites psi[p] and psi[c] with a larger
# common link; the state is unchanged. Returns the number of directions
# added.
function _cbe_expand_bond!(psi, H, p, c, Ep, Ec, maxdim, delta, eps_pre, eps_final)
  lk = commonind(psi[p], psi[c])
  k = dim(lk)
  w = commonind(H[p], H[c])
  d_ex = ceil(Int, (1.0 + delta)*maxdim)
  dtilde = max(0, d_ex - k)
  dpre = max(1, ceil(Int, d_ex/dim(w)))
  dtilde == 0 && return 0
  sp = siteind(psi, p)
  sc = siteind(psi, c)
  op = _cbe_outer_link(psi, p, c)
  oc = _cbe_outer_link(psi, c, p)
  nl = _cbe_dim(op)*dim(sp)
  nr = dim(sc)*_cbe_dim(oc)
  Xt = Ep*psi[p]*H[p]
  X3 = _cbe_array(Xt, [_cbe_prime(op), prime(sp), w, lk], (nl, dim(w), k))
  Yt = psi[c]*H[c]*Ec
  Y3 = _cbe_array(Yt, [lk, w, prime(sc), _cbe_prime(oc)], (k, dim(w), nr))
  Am = _cbe_array(psi[p], [op, sp, lk], (nl, k))
  Cm = _cbe_array(psi[c], [lk, sc, oc], (k, nr))
  At = _cbe_kernel(X3, Y3, Am, Cm, dtilde, dpre, eps_pre, eps_final)
  nt = size(At, 2)
  nt == 0 && return 0
  new = Index(k + nt; tags=tags(lk))
  Pex = hcat(Am, At)
  Cex = vcat(Cm, zeros(eltype(Cm), nt, nr))
  psi[p] = _cbe_itensor(Pex, [op, sp, new])
  psi[c] = _cbe_itensor(Cex, [new, sc, oc])
  return nt
end

_cbe_extend(E, psi, H, i) = E*psi[i]*H[i]*dag(prime(psi[i]))

# Lowest eigenpair of the one-site effective Hamiltonian at site i, from
# psi[i] plus a random vector of relative weight noise (pyitensor's
# _noisy_start), by KrylovKit's Lanczos.
function _cbe_local_solve(L, Ri, H, x0, i, noise, rng)
  if noise > 0
    r = ITensor(randn(rng, Tuple(dim(ix) for ix in inds(x0))...), inds(x0)...)
    x0 = x0 + (noise*norm(x0)/norm(r))*r
  end
  x0 = x0/norm(x0)
  Wi = H[i]
  f(x) = noprime(((L*x)*Wi)*Ri)
  n = dim(x0)
  vals, vecs, _ = _CBE_KK.eigsolve(f, x0, 1, :SR; ishermitian=true,
      tol=1e-13, krylovdim=min(n, 30), maxiter=4, verbosity=0)
  return real(vals[1]), vecs[1]
end

# One left-to-right then right-to-left CBE sweep (pyitensor's _cbe_sweep).
# psi enters with its center on site 1 and R valid for it, and leaves the
# same way. Returns the last local eigenvalue.
function _cbe_sweep!(psi, H, Lenv, Renv, maxdim, cutoff, noise, delta, rng)
  N = length(psi)
  energy = 0.0
  for i in 1:(N - 1)
    _cbe_expand_bond!(psi, H, i + 1, i, Renv[i + 2], Lenv[i - 1], maxdim,
                      delta, 0.0, 0.0)
    Renv[i + 1] = _cbe_extend(Renv[i + 2], psi, H, i + 1)
    energy, theta = _cbe_local_solve(Lenv[i - 1], Renv[i + 1], H, psi[i], i,
                                     noise, rng)
    left = i > 1 ? (commonind(psi[i], psi[i - 1]), siteind(psi, i)) :
                   (siteind(psi, i),)
    U, S, V = ITensors.svd(theta, left...; cutoff=cutoff, maxdim=maxdim,
                  lefttags="Link,l=$i", righttags="Link,l=$i")
    psi[i] = U
    psi[i + 1] = (S*V)*psi[i + 1]
    Lenv[i] = _cbe_extend(Lenv[i - 1], psi, H, i)
  end
  for i in N:-1:2
    _cbe_expand_bond!(psi, H, i - 1, i, Lenv[i - 2], Renv[i + 1], maxdim,
                      delta, 0.0, 0.0)
    Lenv[i - 1] = _cbe_extend(Lenv[i - 2], psi, H, i - 1)
    energy, theta = _cbe_local_solve(Lenv[i - 1], Renv[i + 1], H, psi[i], i,
                                     noise, rng)
    ll = commonind(psi[i], psi[i - 1])
    U, S, V = ITensors.svd(theta, ll; cutoff=cutoff, maxdim=maxdim,
                  lefttags="Link,l=$(i-1)", righttags="Link,l=$(i-1)")
    psi[i] = V
    psi[i - 1] = psi[i - 1]*(U*S)
    Renv[i] = _cbe_extend(Renv[i + 1], psi, H, i)
  end
  return energy
end

# Ground-state DMRG with one-site updates and controlled bond expansion,
# the julia_live counterpart of pyitensor/cbe.py's dmrg_cbe. Every sweep
# runs at maxm and cutoff, as get_gs_dmrg's schedule does; the noise
# goes into the first sweep of a fresh start only, as in dmrg_cbe. Works
# on a copy of psi0 and returns (energy, psi).
function get_gs_cbe(H, psi0; nsweeps=10, cutoff=1e-8, maxm=80, noise=0.0,
                    delta=0.1, fresh=true)
  N = length(psi0)
  N < 2 && error("dmrg_method=\"CBE\" needs at least two sites")
  psi = noprime(copy(psi0), "Link")
  orthogonalize!(psi, 1)
  psi[1] = psi[1]/norm(psi[1])
  rng = Random.MersenneTwister(0)
  Lenv = Dict{Int,ITensor}(0 => ITensor(1.0), -1 => ITensor(1.0))
  Renv = Dict{Int,ITensor}(N + 1 => ITensor(1.0), N + 2 => ITensor(1.0))
  for i in N:-1:2
    Renv[i] = _cbe_extend(Renv[i + 1], psi, H, i)
  end
  energy = 0.0
  for sw in 1:nsweeps
    nz = (fresh && sw == 1) ? float(noise) : 0.0
    energy = _cbe_sweep!(psi, H, Lenv, Renv, maxm, cutoff, nz, delta, rng)
  end
  return energy, psi
end

