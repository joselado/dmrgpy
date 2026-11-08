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
