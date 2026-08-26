"""to_mpo(): turns an AutoMPO (Phase 3) into an actual MPO (Phase 4) --
mirrors mo_terms.h's build_mpo()/toMPO(ampo,{"MaxDim",...,"Exact",false}).

As described in autompo.py's module docstring, this doesn't port ITensor's
own automaton MPO-compression algorithm: each HTerm becomes its own exact,
trivial (bond dimension 1) MPO, and all of them are combined via
mpsalgebra.sum_many() -- one exact, K-way block-diagonal concatenation
(pure array placement, no SVD) followed by exactly the bidirectional
compression pass described next.

Concatenation alone (with no compression at all) is *not* enough to reach
the true minimal bond dimension, though -- confirmed directly, not just
assumed: a plain nearest-neighbor Heisenberg chain's MPO came out with
bond dimension equal to its term count (39 at N=14) instead of the
well-known constant ~5. A one-directional (left-to-right only) truncating
sweep doesn't fix this either: it only ever compresses relative to the
*bonds already finalized to its left*, and can't see that (for instance)
the "still need to place every remaining identity" tail is the same
redundant structure repeated at every term. A bidirectional compression
pass (position() rightward then leftward, both truncating) after
concatenation fixes this completely -- confirmed directly on the same
case, 39 -> 5 -- because sweeping both directions lets SVD see the
*global* redundancy, not just what's accumulated so far in one direction.
Doing that once, after a single exact concatenation of every term, keeps
the SVD work bounded to O(1) sweeps regardless of term count T, instead of
sum_many()'s pairwise predecessor (repeated 2-way sum()), which paid a
full truncating SVD sweep on every one of T-1 folds -- see sum_many()'s
own docstring (mpsalgebra.py) for the T-1-wasted-sweeps argument and a
benchmark comparing the two.
"""

import numpy as np

from . import backend as bk

from .index import Index
from .mpscontainer import MPO
from .mpsalgebra import sum_many as _mps_sum_many
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
        tensors.append(ITensor(tuple(inds), bk.zeros(shape)))
        prev_link = right_link
    mpo = MPO(tensors)
    mpo.center = 1
    return mpo


def to_mpo(ampo, cutoff=0.0, maxdim=None):
    if not ampo.terms:
        return _zero_mpo(ampo.sites)
    term_mpos = [_term_to_mpo(term, ampo.sites) for term in ampo.terms]
    result = _mps_sum_many(term_mpos, cutoff=cutoff, maxdim=maxdim)
    if result.length() > 1:
        # One-directional per-pairwise-sum compression alone leaves real
        # redundancy on the table (see this module's docstring) -- a final
        # there-and-back sweep reaches the true minimal bond dimension.
        result.position(result.length(), cutoff=cutoff, maxdim=maxdim)
        result.position(1, cutoff=cutoff, maxdim=maxdim)
    return result
