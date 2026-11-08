using ITensors
using ITensorMPS
using LinearAlgebra

# Four-point correlator tensor <psi| O1_i O2_j O3_k O4_l |psi> of fermionic
# operators on a chain of spinless Fermion sites, with environment reuse.
#
# In the Jordan-Wigner representation every factor is a tensor product,
# O_f at site s_f = F_1 ... F_{s_f - 1} a_f(s_f), so the product of the
# four factors is again a tensor product over sites, and the local matrix
# at site m is the ordered product over f = 1..4 of
#   identity     if s_f < m  (the factor is already placed),
#   a_f          if s_f == m (the factor sits here),
#   F            if s_f > m  (its string still runs through m).
# No reordering sign is needed: the order of the factors is kept site by
# site. Repeated indices (i == j, i == k, ...) come out of the same rule.
#
# The tensor is filled by a depth-first walk over the sites: the state of
# the walk at site m is the set of factors placed so far and their sites,
# so every tuple that shares a prefix shares its left environment. A tuple
# closes on the site where its last factor is placed, against the
# identity right environment of the remaining sites. The left and right
# identity environments are built once, so the result is the raw
# <psi|O|psi> for any MPS, normalized or not, in any gauge.

# local matrix M[s',s] of an operator name on site index idx
function _fp_matrix(name::String, idx::Index)
  return Matrix{ComplexF64}(array(op(name, idx), prime(idx), idx))
end

# local operator at site m for the per-factor states in codes
# (0 placed before, 1 here, 2 later); nothing means the identity
function _fp_local_op!(cache, sites, fnames, m, codes)
  key = (m, codes)
  haskey(cache, key) && return cache[key]
  idx = sites[m]
  d = dim(idx)
  M = Matrix{ComplexF64}(I, d, d)
  isid = true
  for f in 1:length(codes)
    c = codes[f]
    if c == 1
      M = M * _fp_matrix(fnames[f], idx)
      isid = false
    elseif c == 2
      M = M * _fp_matrix("F", idx)
      isid = false
    end
  end
  # an even number of F strings and no operator is the identity
  if isid || M ≈ Matrix{ComplexF64}(I, d, d)
    cache[key] = nothing
  else
    cache[key] = ITensor(M, prime(idx), dag(idx))
  end
  return cache[key]
end

function _fp_fold(E, psi, bra, braL, O, m)
  if O === nothing
    return E * psi[m] * braL[m]
  end
  return E * psi[m] * O * bra[m]
end

function _fp_walk!(out, E, m, pos, U, n, psi, bra, braL, R, cache, sites, fnames)
  nf = length(pos)
  # every subset S of the unplaced factors U can be placed on site m
  S = U
  while true
    if m < n || S == U # the last site must place every factor left
      codes = ntuple(f -> pos[f] > 0 ? 0 : (((S >> (f-1)) & 1) == 1 ? 1 : 2), nf)
      O = _fp_local_op!(cache, sites, fnames, m, codes)
      T = _fp_fold(E, psi, bra, braL, O, m)
      for f in 1:nf
        ((S >> (f-1)) & 1) == 1 && (pos[f] = m)
      end
      if S == U # every factor placed: close against the right environment
        out[pos...] += scalar(T * R[m+1])
      else
        _fp_walk!(out, T, m+1, pos, U & ~S, n, psi, bra, braL, R, cache,
                  sites, fnames)
      end
      for f in 1:nf
        ((S >> (f-1)) & 1) == 1 && (pos[f] = 0)
      end
    end
    S == 0 && break
    S = (S - 1) & U
  end
  return nothing
end

"""four_point_sweep(psi, fnames) -> Array{ComplexF64,4}

out[i,j,k,l] = <psi| O1_i O2_j O3_k O4_l |psi> for the operator names
fnames (each a fermionic operator, e.g. ["Cdag","C","Cdag","C"]) on a
chain of Fermion sites, Jordan-Wigner strings included."""
function four_point_sweep(psi::MPS, fnames::Vector{String})
  n = length(psi)
  nf = length(fnames)
  sites = siteinds(psi)
  bra = [dag(prime(psi[m])) for m in 1:n]
  braL = [dag(prime(psi[m], "Link")) for m in 1:n]
  L = Vector{ITensor}(undef, n+1) # L[m]: identity environment of sites 1..m-1
  L[1] = ITensor(1.0)
  for m in 1:n
    L[m+1] = L[m] * psi[m] * braL[m]
  end
  R = Vector{ITensor}(undef, n+1) # R[m]: identity environment of sites m..n
  R[n+1] = ITensor(1.0)
  for m in n:-1:1
    R[m] = psi[m] * braL[m] * R[m+1]
  end
  out = zeros(ComplexF64, ntuple(_ -> n, nf)...)
  cache = Dict{Any,Any}()
  pos = zeros(Int, nf)
  full = (1 << nf) - 1
  for mn in 1:n # site of the first placed factor
    S = full
    while S != 0
      if mn == n && S != full # the last site must place every factor
        S = (S - 1) & full
        continue
      end
      codes = ntuple(f -> ((S >> (f-1)) & 1) == 1 ? 1 : 2, nf)
      O = _fp_local_op!(cache, sites, fnames, mn, codes)
      T = _fp_fold(L[mn], psi, bra, braL, O, mn)
      for f in 1:nf
        ((S >> (f-1)) & 1) == 1 && (pos[f] = mn)
      end
      if S == full
        out[pos...] += scalar(T * R[mn+1])
      else
        _fp_walk!(out, T, mn+1, pos, full & ~S, n, psi, bra, braL, R, cache,
                  sites, fnames)
      end
      fill!(pos, 0)
      S = (S - 1) & full
    end
  end
  return out
end
