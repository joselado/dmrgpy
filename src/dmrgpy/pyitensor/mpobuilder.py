"""to_mpo(): turns an AutoMPO (Phase 3) into an actual MPO (Phase 4) --
mirrors mo_terms.h's build_mpo()/toMPO(ampo,{"MaxDim",...,"Exact",false}).

As described in autompo.py's module docstring, this doesn't port ITensor's
own automaton MPO-compression algorithm: each HTerm becomes its own exact,
trivial (bond dimension 1) MPO, and all of them are summed together via
mpsalgebra.sum(), which SVD-compresses after every pairwise sum.

That per-step compression alone is *not* enough to reach the true minimal
bond dimension, though -- confirmed directly, not just assumed: a plain
nearest-neighbor Heisenberg chain's MPO came out with bond dimension
growing linearly in N (39 at N=14) instead of the well-known constant ~5,
because each one-directional (left-to-right) SVD sweep only ever
compresses relative to the *bonds already finalized to its left*, and
can't see that (for instance) the "still need to place every remaining
identity" tail is the same redundant structure repeated at every one of
those intermediate sums. A single extra bidirectional compression pass
(position() rightward then leftward, both truncating) after all terms are
summed fixes this completely -- confirmed directly on the same case,
39 -> 5 -- because sweeping both directions lets SVD see the *global*
redundancy, not just what's accumulated so far in one direction. Doing
that once at the end (rather than compressing bidirectionally after every
single pairwise sum) keeps intermediate work bounded without sacrificing
the final result's bond dimension.
"""

import numpy as np

from .index import Index
from .mpscontainer import MPO
from .mpsalgebra import sum as _mps_sum
from .tensor import ITensor


def _term_to_mpo(term, sites):
    n = sites.length()
    mats = term.resolve(sites)  # standard-convention (dim,dim) matrices, index 0 = site 1
    tensors = []
    prev_link = None
    for i in range(1, n + 1):
        s = sites.si(i)
        stored = mats[i - 1].T  # std (out,in) -> this engine's (in,out) storage convention
        left_link = prev_link
        right_link = Index(1, tags="Link,l={}".format(i)) if i < n else None

        inds = ([left_link] if left_link else []) + [s, s.prime(1)] + ([right_link] if right_link else [])
        shape = [1] if left_link else []
        shape += [s.dim, s.dim]
        if right_link:
            shape += [1]
        arr = stored.reshape(tuple(shape))
        tensors.append(ITensor(tuple(inds), arr))
        prev_link = right_link

    tensors[0] = tensors[0] * term.coef
    mpo = MPO(tensors)
    mpo.center = 1
    return mpo


def _zero_mpo(sites):
    """The zero operator, as a trivial bond-dimension-1 MPO with an
    all-zero matrix at every site -- the mathematically sensible reading
    of "a sum of zero terms", needed because dmrgpy's own backend-agnostic
    code (e.g. algebra/arnolditk.py's Arnoldi orthogonalize(), which
    multiplies an MPS by `coefficient*multioperator.identity()`) can
    legitimately produce an empty AutoMPO whenever that coefficient gets
    filtered to (numerically) zero -- confirmed directly: a 2-site
    Spinful_Fermionic_Chain's very first Arnoldi orthogonalization step
    hit exactly this, since there's nothing yet to project out."""
    n = sites.length()
    tensors = []
    prev_link = None
    for i in range(1, n + 1):
        s = sites.si(i)
        left_link = prev_link
        right_link = Index(1, tags="Link,l={}".format(i)) if i < n else None
        inds = ([left_link] if left_link else []) + [s, s.prime(1)] + ([right_link] if right_link else [])
        shape = tuple(ind.dim for ind in inds)
        tensors.append(ITensor(tuple(inds), np.zeros(shape, dtype=complex)))
        prev_link = right_link
    mpo = MPO(tensors)
    mpo.center = 1
    return mpo


def to_mpo(ampo, cutoff=0.0, maxdim=None):
    if not ampo.terms:
        return _zero_mpo(ampo.sites)
    result = None
    for term in ampo.terms:
        term_mpo = _term_to_mpo(term, ampo.sites)
        result = term_mpo if result is None else _mps_sum(result, term_mpo, cutoff=cutoff, maxdim=maxdim)
    if result.length() > 1:
        # One-directional per-pairwise-sum compression alone leaves real
        # redundancy on the table (see this module's docstring) -- a final
        # there-and-back sweep reaches the true minimal bond dimension.
        result.position(result.length(), cutoff=cutoff, maxdim=maxdim)
        result.position(1, cutoff=cutoff, maxdim=maxdim)
    return result
