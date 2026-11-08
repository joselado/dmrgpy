using ITensors
using ITensorMPS

# Ground-state and expectation-value helpers of the julia_live backend
# (mpsjulialive/vev.py, mpsjulialive/groundstate.py, entropytk/
# correlationentropy.py). Loaded last by juliasession.py's initialize().


# Copy of an MPS without the stale prime level that operator algebra leaves
# on its Link indices (see groundstate.py's get_gs_dmrg note), so that the
# bra and the MPO below meet on clean indices.
function gs_extra_clean(psi)
  return noprime(copy(psi),"Link")
end


# <a|M|b>, contracted exactly: no MPO-MPS application, so no truncation to
# maxdim or cutoff. vev() ends with it.
function gs_extra_aMb(a,M,b)
  return inner(gs_extra_clean(a)',M,gs_extra_clean(b))
end


# M|psi>, truncated to maxdim and to the relative cutoff, with plain site
# indices on the result.
function gs_extra_apply(M,psi,maxdim,cutoff)
  return apply(M,gs_extra_clean(psi);maxdim=maxdim,cutoff=cutoff,
               alg="densitymatrix")
end


# <Cdag_i C_j Cdag_k C_l> for every (i,j,k,l), one exact MPO per tuple, the
# counterpart of pyitensor/chain.py's four_correlation_tensor. With
# accelerate, a tuple after its conjugate partner (l,k,j,i) in lexicographic
# order is not computed, and is filled with the complex conjugate. Indices
# are 0-based in the returned array as seen from Python (1-based here).
function four_correlation_tensor(sites,psi,accelerate)
  n = length(sites)
  out = zeros(ComplexF64,n,n,n,n)
  phi = gs_extra_clean(psi)
  for i in 1:n, j in 1:n, k in 1:n, l in 1:n
    current = (i,j,k,l)
    partner = (l,k,j,i)
    if accelerate && current > partner
      continue
    end
    os = OpSum()
    os += 1.0,"Cdag",i,"C",j,"Cdag",k,"C",l
    M = MPO(os,sites)
    c = inner(phi',M,phi)
    out[i,j,k,l] = c
    if accelerate && current != partner
      out[l,k,j,i] = conj(c)
    end
  end
  return out
end
