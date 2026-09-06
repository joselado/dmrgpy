"""to_mpo(): turns an AutoMPO (Phase 3) into an actual MPO (Phase 4) --
mirrors mo_terms.h's build_mpo()/toMPO(ampo,{"MaxDim",...,"Exact",false}).

How the MPO is built, and why it changed
----------------------------------------
The operator is assembled directly as a **finite-state machine** over the
partial products of the terms (`_automaton_mpo` below), the same shape
ITensor's own `toMPO(...,{"Exact",true})` produces, and then handed one
bidirectional truncating sweep so the caller's `cutoff`/`maxdim` are still
honoured exactly as before.

This replaced the original construction, which is worth recording because
its failure mode was counter-intuitive. That version built one exact,
trivial (bond dimension 1) MPO per HTerm and combined all T of them with
`mpsalgebra.sum_many()` -- one exact, K-way block-diagonal concatenation
(pure array placement, no SVD) followed by the same bidirectional
compression pass. The *final* bond dimension it reached was right (39 -> 5
on a nearest-neighbour Heisenberg chain at N=14, the well-known constant),
and the compression pass genuinely was needed to get there: concatenation
alone left bond dimension equal to the term count, and a one-directional
sweep didn't fix it either, because it only ever compresses relative to
the bonds already finalized to its left and cannot see that the "still
need to place every remaining identity" tail is the same redundant
structure repeated at every term. Sweeping both directions lets SVD see
the *global* redundancy. All of that was confirmed directly and all of it
is still true.

**The cost was the intermediate, not the answer.** The concatenated MPO
carries bond dimension ~T = 3(L-1) for a nearest-neighbour chain *before*
compression, so the compression sweep ran O(L) truncating SVDs on matrices
of size O(L): O(L^4) overall. Measured on an S=1/2 Heisenberg chain
(`maxm=30`, `nsweeps=15`, threads pinned), one `gs_energy(mode="DMRG")`
took 1.08 s at L=20, 8.29 s at L=40 and 239 s at L=100 -- i.e. ~L^3.7
where fixed maxm and fixed nsweeps should give ~L -- and at L=100 `to_mpo`
was 95% of the whole "ground state calculation", 76% of it inside
`svd._svd_truncated`. The DMRG itself was ~13 s and scaled linearly, i.e.
correctly. See docs/pip_install_and_pyitensor_performance_plan.md for the
full profile.

The automaton reaches the same minimal bond dimension without ever
building the O(L)-dimensional intermediate. The machine itself is O(L):
one bond at a time, with the compression sweep now running at the
*final*, O(1), bond dimension. Be precise about the total, though -- it is
O(L^2), not O(L), because `HTerm.resolve()` spells every term out over
all L sites and there are O(L) terms. That term is cheap (small dense
per-site matrices, ~30k of them at L=100) and is nowhere near dominant at
the sizes measured, but it is the next wall if one ever appears; the
O(L^4) it replaced is gone either way.

The machine
-----------
Each term, once `HTerm.resolve()` has spelled it out as one matrix per
site (Jordan-Wigner F-strings included), occupies sites `[first,last]` --
the first and last site where that matrix is not the identity. Across a
bond it is therefore in one of three situations, which are exactly the
machine's states:

* **I** -- not started yet (identity everywhere to the left);
* **F** -- finished (identity everywhere to the right);
* one **partial** state per *distinct* left-partial product among the
  terms straddling this bond.

Sharing partial states between terms is what compresses: two terms with
the same matrices on sites 1..k share their state at bond k and branch
afterwards, which is why the nearest-neighbour Heisenberg chain comes out
at 2 + 3 = 5 and not at its term count. It creates no spurious paths --
sharing a state *means* the prefixes are identical, so following one
term's prefix into another's suffix just reproduces that other term.

Two rules make the bookkeeping correct rather than merely plausible:

* a term's **coefficient goes on its transition into F**, never earlier,
  so terms differing only by a coefficient still share every partial
  state;
* transitions into F **accumulate** (`+=`) while structural transitions
  are **assigned**. Two syntactically identical terms trace the very same
  path, so their coefficients must sum -- but their shared structural
  transitions must not be written twice. (The two kinds never target the
  same matrix entry: F is a column no structural transition writes, and
  no term's transition ever leaves the F row.)
"""

import numpy as np

from . import backend as bk

from .index import Index
from .mpscontainer import MPO
from .mpsalgebra import sum_many as _mps_sum_many
from .tensor import ITensor


def _term_to_mpo(term, sites):
    """One HTerm as an exact bond-dimension-1 MPO. No longer on the
    production path (see this module's docstring) -- kept because
    `_sum_of_term_mpos` below, the construction this module used to use,
    is an independent reference implementation to check the automaton
    against (tests/test_mpo_automaton_builder.py)."""
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


def _sum_of_term_mpos(ampo, cutoff=0.0, maxdim=None):
    """The pre-automaton construction, kept as a reference implementation:
    one bond-dimension-1 MPO per term, block-diagonally concatenated. Its
    intermediate bond dimension is the term count, which is what made it
    O(L^4); `_automaton_mpo` is the production path now."""
    term_mpos = [_term_to_mpo(term, ampo.sites) for term in ampo.terms]
    return _mps_sum_many(term_mpos, cutoff=cutoff, maxdim=maxdim)


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


def _is_identity(mat):
    """Exact (bitwise) test, not a tolerance: every untouched site's matrix
    comes straight from `matrix("Id")`/a non-fermionic site's trivial "F",
    both of which are literally np.eye. A matrix that is only *numerically*
    the identity is simply treated as part of the term's support, which
    costs at most one extra automaton state and is then removed by the
    compression sweep -- conservative, never wrong."""
    return np.array_equal(mat, np.eye(mat.shape[0]))


def _term_spans(resolved):
    """(first,last) site of each term's non-identity support, 1-based.

    A term that resolves to a multiple of the identity everywhere (a bare
    coefficient, which arnolditk.py's orthogonalization really does
    produce) has no support at all; it is reported as (1,1) so it enters
    the machine as a single I->F transition at site 1, contributing
    coef*Id to the whole chain, which is exactly what it means.

    Note this deliberately reads the *resolved* matrices rather than
    `term.ops`: a term with an odd number of fermionic factors carries its
    Jordan-Wigner string all the way to the right end of the chain, so its
    support genuinely runs past its last named operator."""
    spans = []
    for mats in resolved:
        nz = [i for i, m in enumerate(mats) if not _is_identity(m)]
        spans.append((nz[0] + 1, nz[-1] + 1) if nz else (1, 1))
    return spans


def _automaton_mpo(ampo):
    """Build the MPO directly as the finite-state machine described in this
    module's docstring. Exact: no truncation happens here at all (the
    caller's cutoff/maxdim are applied afterwards, by to_mpo's sweep)."""
    sites = ampo.sites
    n = sites.length()
    terms = ampo.terms
    resolved = [term.resolve(sites) for term in terms]
    spans = _term_spans(resolved)

    if n == 1:
        # No bonds to run a machine over: the MPO is the one-site sum.
        s = sites.si(1)
        arr = np.zeros((s.dim, s.dim), dtype=complex)
        for term, mats in zip(terms, resolved):
            arr += term.coef * mats[0].T
        mpo = MPO((ITensor((s, s.prime(1)), arr),))
        mpo.center = 1
        return mpo

    state = [0] * len(terms)  # each term's row index at the previous bond
    row_dim = 1               # bond 0 carries the I state alone
    tensors = []
    prev_link = None
    for k in range(1, n + 1):
        s = sites.si(k)
        d = s.dim
        eye = np.eye(d)
        last_site = (k == n)
        # Column layout at bond k: I at 0 and F at 1, partials from 2 on --
        # except at the right edge, where only F survives.
        col_F = 0 if last_site else 1
        partial_cols = {}
        writes = []  # (row, col, standard-convention matrix, accumulate?)
        if not last_site:
            writes.append((0, 0, eye, False))        # I -> I: nothing started yet
        if k >= 2:
            writes.append((1, col_F, eye, False))    # F -> F: everything already placed
        for t, term in enumerate(terms):
            first, last = spans[t]
            if k < first or k > last:
                continue  # this term sits in I (or in F) and rides the channels above
            row = 0 if k == first else state[t]
            mat = resolved[t][k - 1]
            if k == last:
                writes.append((row, col_F, term.coef * mat, True))
            else:
                # The state is "which partial product got us here", so
                # the key is the incoming state plus this site's matrix.
                # Normalized to a contiguous complex buffer first: raw
                # .tobytes() is stride- and dtype-sensitive, so two equal
                # matrices stored differently would silently fail to share
                # a state (harmless -- a bigger MPO the sweep then
                # compresses -- but silently losing the whole point).
                key = (row, np.ascontiguousarray(mat, dtype=complex).tobytes())
                col = partial_cols.get(key)
                if col is None:
                    col = 2 + len(partial_cols)
                    partial_cols[key] = col
                writes.append((row, col, mat, False))
                state[t] = col
        col_dim = 1 if last_site else 2 + len(partial_cols)

        W = np.zeros((row_dim, d, d, col_dim), dtype=complex)
        for row, col, mat, accumulate in writes:
            # std (out,in) -> this engine's (in,out) storage convention,
            # exactly as _term_to_mpo does it
            if accumulate:
                W[row, :, :, col] += mat.T
            else:
                W[row, :, :, col] = mat.T

        left_link = prev_link
        right_link = None if last_site else Index(col_dim, tags="Link,l={}".format(k))
        inds = ([left_link] if left_link is not None else []) + [s, s.prime(1)] \
            + ([right_link] if right_link is not None else [])
        shape = ([row_dim] if left_link is not None else []) + [d, d] \
            + ([col_dim] if right_link is not None else [])
        tensors.append(ITensor(tuple(inds), W.reshape(tuple(shape))))
        prev_link = right_link
        row_dim = col_dim

    mpo = MPO(tensors)
    # Not canonical yet, and position() does not need it to be: an SVD
    # split is exact regardless of the tensor's prior gauge, so the sweep
    # to_mpo runs next canonicalizes as it goes (see _shift_right's
    # docstring in mpscontainer.py). _term_to_mpo has always done the same.
    mpo.center = 1
    return mpo


def to_mpo(ampo, cutoff=0.0, maxdim=None):
    if not ampo.terms:
        return _zero_mpo(ampo.sites)
    result = _automaton_mpo(ampo)
    if result.length() > 1:
        # The machine is exact but not canonical and not truncated: this is
        # where the caller's cutoff/maxdim actually get applied, and where
        # any redundancy the partial-prefix sharing cannot see (suffix
        # sharing, linearly dependent channels) is squeezed out. Cheap now
        # that the incoming bond dimension is O(1) rather than O(T).
        result.position(result.length(), cutoff=cutoff, maxdim=maxdim)
        result.position(1, cutoff=cutoff, maxdim=maxdim)
    return result
