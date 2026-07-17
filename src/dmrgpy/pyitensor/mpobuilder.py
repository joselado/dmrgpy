"""to_mpo(): turns an AutoMPO (Phase 3) into an actual MPO (Phase 4) --
mirrors mo_terms.h's build_mpo()/toMPO(ampo,{"MaxDim",...,"Exact",false}).

As described in autompo.py's module docstring, this doesn't port ITensor's
own automaton MPO-compression algorithm: each HTerm becomes its own exact,
trivial (bond dimension 1) MPO, and all of them are summed together via
mpsalgebra.sum(), which does the actual SVD compression down to
cutoff/maxdim. Simpler to get right than a state-merging automaton, and
asymptotically fine for the term counts dmrgpy's own Hamiltonians have.
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
    return result
