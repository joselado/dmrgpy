"""svd(): the one tensor-factorization primitive the rest of this engine is
built on -- every truncation in the library (MPS canonicalization, sum(),
applyMPO(), nmultMPO(), the DMRG/TDVP local update) is "reshape into a
matrix across some index partition, SVD, keep the largest singular values
up to Cutoff/MaxDim, reshape back". mpscpp3/chain_session.h only calls it
directly once (bond_entropy()); every other MPS/MPO-level truncation goes
through this same function via mps.py/mpo.py in a later phase.
"""

import numpy as np

from .index import Index
from .tensor import ITensor, _find


class Spectrum:
    """The truncated singular-value spectrum of one svd() call. `probs` are
    the kept singular values squared and normalized against the *original*
    (pre-truncation) total, so sum(probs) + truncerr == 1 -- mirrors
    ITensor's Spectrum::eigs()/truncerr(), which bond_entropy() uses
    directly to compute the von Neumann entropy of a bond."""

    def __init__(self, singular_values, probs, truncerr):
        self.singular_values = singular_values
        self._probs = probs
        self.truncerr = truncerr

    def eigs(self):
        return self._probs


def _truncate(s, cutoff, maxdim, mindim):
    """s: singular values, sorted descending. Returns (keep, discarded_weight)
    -- how many of the largest to keep, and the normalized weight (sum of
    p_i over the dropped tail) that was discarded. Mirrors ITensor's own
    truncation rule: drop the smallest singular values first, stopping as
    soon as either mindim is reached, or dropping the next one would exceed
    `cutoff` of discarded weight -- except maxdim is a hard cap that's
    enforced regardless of cutoff."""
    n = len(s)
    p = s.astype(float) ** 2
    total = p.sum()
    if total <= 0:
        return max(1, mindim), 0.0
    p = p / total
    keep = n
    discarded = 0.0
    while keep > mindim:
        over_maxdim = maxdim is not None and keep > maxdim
        if not over_maxdim and discarded + p[keep - 1] > cutoff:
            break
        discarded += p[keep - 1]
        keep -= 1
    return keep, discarded


def svd(T, left_inds, cutoff=0.0, maxdim=None, mindim=1, tags="Link"):
    """Split ITensor T into U * S * V (contracting U*S*V reconstructs T up
    to the requested truncation), grouping `left_inds` onto U and every
    other index of T onto V. Returns (U, S, V, spectrum) -- mirrors
    ITensor's own svd(T,U,S,V) (there an in/out-param triple plus a
    returned Spectrum; here a 4-tuple since Python has no out-params).

    U and V each get their own freshly minted bond index (both dimension
    `keep`, both tagged `tags`) rather than sharing one -- so U's columns
    and V's rows are each independently orthonormal, and callers can absorb
    S into either side (`U*S` for a right-canonical partner, `S*V` for a
    left-canonical one) depending on which way an MPS sweep is moving,
    exactly as ITensor's own three-tensor split allows.
    """
    left_inds = list(left_inds)
    for ind in left_inds:
        if not T.hasindex(ind):
            raise ValueError("svd: {} is not an index of {}".format(ind, T))
    right_inds = [ind for ind in T.inds if _find(left_inds, ind) is None]

    order = left_inds + right_inds
    arr = T.transpose_to(order)
    ldim = int(np.prod([ind.dim for ind in left_inds], dtype=int)) if left_inds else 1
    rdim = int(np.prod([ind.dim for ind in right_inds], dtype=int)) if right_inds else 1
    mat = arr.reshape(ldim, rdim)

    U, S, Vh = np.linalg.svd(mat, full_matrices=False)
    keep, discarded = _truncate(S, cutoff, maxdim if maxdim else None, max(1, mindim))

    probs_full = (S.astype(float) ** 2)
    total = probs_full.sum()
    probs = probs_full[:keep] / total if total > 0 else probs_full[:keep]

    bond_u = Index(keep, tags=tags)
    bond_v = Index(keep, tags=tags)
    left_shape = tuple(ind.dim for ind in left_inds) + (keep,)
    right_shape = (keep,) + tuple(ind.dim for ind in right_inds)

    Utensor = ITensor(tuple(left_inds) + (bond_u,), U[:, :keep].reshape(left_shape))
    Vtensor = ITensor((bond_v,) + tuple(right_inds), Vh[:keep, :].reshape(right_shape))
    Stensor = ITensor((bond_u, bond_v), np.diag(S[:keep].astype(complex)))
    spectrum = Spectrum(S[:keep], probs, discarded)
    return Utensor, Stensor, Vtensor, spectrum
