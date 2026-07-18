# Dense-vector/matrix <-> MPS/MPO conversion helpers, for testing only (not
# part of the pyitensor package itself -- chain_session.h has no equivalent
# "to dense" API, so there's no production need for this outside tests).
import numpy as np

from dmrgpy.pyitensor import MPS, MPO, ITensor, Index, svd


def mps_to_dense(mps):
    """Contract an MPS to a dense vector, site 1 most significant -- the
    same convention as edtk/one2many.py's own np.kron loop and
    autompo.py's dense_matrix()."""
    n = mps.length()
    phys = [next(ind for ind in mps.A(i).inds if ind.hastags("Site")) for i in range(1, n + 1)]
    vec = mps.A(1)
    for i in range(2, n + 1):
        vec = vec * mps.A(i)
    arr = vec.transpose_to(phys)
    return arr.reshape(-1)


def mpo_to_dense(mpo):
    """Contract an MPO to a dense (dim,dim) standard-convention matrix."""
    n = mpo.length()
    phys_in = [next(ind for ind in mpo.A(i).inds if ind.hastags("Site") and ind.plev == 0)
               for i in range(1, n + 1)]
    phys_out = [next(ind for ind in mpo.A(i).inds if ind.hastags("Site") and ind.plev == 1)
                for i in range(1, n + 1)]
    T = mpo.A(1)
    for i in range(2, n + 1):
        T = T * mpo.A(i)
    arr = T.transpose_to(phys_in + phys_out)
    din = int(np.prod([p.dim for p in phys_in])) if phys_in else 1
    dout = int(np.prod([p.dim for p in phys_out])) if phys_out else 1
    arr = arr.reshape(din, dout)
    return arr.T  # (in,out) storage -> standard (out,in) convention


def dense_to_mps(vector, sites, cutoff=0.0, maxdim=None):
    """Exact (up to cutoff/maxdim) MPS for a given dense state vector, via
    sequential SVD splitting from the left -- an independent construction
    path from randomMPS()/DMRG, useful as a known-correct reference."""
    n = sites.length()
    dims = [sites.dim(i) for i in range(1, n + 1)]
    phys_inds = [sites.si(i) for i in range(1, n + 1)]
    arr = np.asarray(vector, dtype=complex).reshape(dims)
    cur = ITensor(tuple(phys_inds), arr)

    tensors = []
    left_link = None
    for i in range(1, n):
        left_inds = ([left_link] if left_link else []) + [phys_inds[i - 1]]
        U, S, V, spec = svd(cur, left_inds, cutoff=cutoff, maxdim=maxdim)
        tensors.append(U)
        cur = S * V
        left_link = U.inds[-1]
    tensors.append(cur)

    mps = MPS(tensors)
    mps.center = n
    return mps
