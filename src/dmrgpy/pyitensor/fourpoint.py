"""The four-point tensor <Cdag_i C_j Cdag_k C_l> as batched GEMMs.

Why this exists
---------------
`Chain.four_correlation_tensor_sweep` and `Chain.four_correlation_tensor_fold`
both evaluate one tuple's worth of MPS transfer at a time: a chain of
`chi x chi` by `chi x d x chi` contractions, several hundred thousand of
them at n=30, each far too small to keep a BLAS call -- let alone a GPU --
busy. Measured on a 30-site spinless chain at maxm=20, that is 72 s on the
compiled ITensor v3 backend and ~190 s on `itensor_version="python"`, and
in both cases the *repeated-index* tuples (the O(n^3) of them whose four
indices are not pairwise distinct) cost more than the O(n^4) distinct ones,
because the sweep shares environments across tuples and the repeated-index
fold does not.

This module computes the same tensor with the same arithmetic, reorganized
so that every contraction is one large GEMM over a *batch* of environments
that differ only in which sites their operators sit on. Two things fall out
of that reorganization, and they are the whole point:

* environments are shared across every tuple that agrees on a prefix --
  which subsumes the sweep's own environment reuse *and* extends it to the
  repeated-index tuples the sweep left as per-tuple folds;
* the batch dimension is what makes the calculation worth putting on a
  device at all. `docs/pyitensor_gpu_port_plan.md` Sec. 8 measured the
  eager dispatch floor at ~0.35 ms per array operation, so 164430 separate
  chi=20 contractions would cost ~57 s of pure dispatch before any
  arithmetic. Batched, the same work is ~n GEMMs per trie node.

The trie
--------
Write a tuple's operator content in site-sorted order. It becomes a
sequence of at most four *ranks*: rank r is the r-th smallest distinct site
the four factors occupy, carrying the product of whichever factors sit
there. Two tuples that agree on their first r (matrix, parity) steps have
the same partial environment, whatever sites those steps landed on -- so
the set of all n_modes^4 tuples collapses onto a trie of a few dozen
*matrix* sequences, and each trie node holds one array of environments,
batched over the site combinations that realize it.

A node's batch is therefore `(B, chi, chi)` with `B` up to C(n, r), and one
site of the sweep is, per node, exactly two GEMMs:

    X[b,j,p,k]  = sum_i  E[b,i,j] A_op[i,p,k]
    E'[b,k,l]   = sum_jp X[b,j,p,k] conj(A)[j,p,l]

both reshaped into a single `(B*chi, chi) x (chi, d*chi)` matmul, which is
what NumPy's BLAS and XLA both want.

Conventions (deliberately inherited, not re-derived)
----------------------------------------------------
Every convention here is `Chain.four_correlation_tensor_fold`'s, because
that is the one of the two implementations that already covers *all*
tuples with a single rule -- distinct and repeated alike -- and is
validated against the AutoMPO+`to_mpo`+`inner` pipeline:

* operator matrices are `site_type.matrix(name)`, indexed `[ket, bra]`,
  applied as `tensordot(A, mat, ([1],[0])).transpose(0,2,1)`;
* factors sharing a site compose in increasing factor-position order, so
  the leftmost factor of `Cdag_i C_j Cdag_k C_l` is applied first;
* the Jordan-Wigner string reduces to one rule covering gap sites and
  operator sites alike: a site carries an extra `F` iff the number of
  factors placed up to and including it is odd (`fold()`'s `carry != odd`
  test, which is that statement written incrementally);
* the cross-site reordering sign is `_four_pt_site_sort_sign` of the
  tuple's own site list, which depends only on the rank assignment and is
  therefore constant across every site combination realizing one shape.

The one thing that is *not* inherited is the gauge. The fold re-runs
`psi.position(mn)` per group so it can start from an identity left
environment; a single left-to-right batched pass cannot re-gauge per tuple,
so this uses `psi.position(1)` once and grows the operator-free left
environment `E0` alongside the trie. Closing is still a plain trace,
because `position(1)` leaves every site to the right right-canonical.
"""

import itertools
from functools import lru_cache
from math import comb

import numpy as np

from . import backend as bk
from .tensor import commonIndex
from .sites.base import is_fermionic


def _arrays_lpr(psi):
    """Every site tensor as a `(chi_left, d, chi_right)` array in the
    active backend's namespace.

    `chain._mps_arrays_lpr` does the same thing but finishes with
    `np.ascontiguousarray`, which on a device backend is a silent
    device->host copy per site -- the exact round trip
    `docs/pyitensor_gpu_port_plan.md` Sec. 5 exists to avoid. Here the
    arrays feed GEMMs that must stay where the MPS lives, so the copy is
    replaced by the backend's own `asarray`."""
    xp = bk.xp()
    n = psi.length()
    out = []
    for i in range(1, n + 1):
        T = psi.A(i)
        phys = next(ind for ind in T.inds if ind.hastags("Site"))
        left = commonIndex(psi.A(i - 1), T) if i > 1 else None
        right = commonIndex(T, psi.A(i + 1)) if i < n else None
        order = [ind for ind in (left, phys, right) if ind is not None]
        arr = T.transpose_to(order)
        if left is None:
            arr = arr.reshape((1,) + arr.shape)
        if right is None:
            arr = arr.reshape(arr.shape + (1,))
        out.append(xp.asarray(arr))
    return out


def _transfer(E, A_op, A_conj):
    """One site of the transfer, batched over `E`'s leading axis.

    `E` is `(B, ket, bra)`, `A_op` the ket tensor with this site's operator
    already folded in, `A_conj` the plain conjugated ket tensor. Written as
    two reshaped matmuls rather than `einsum`: at these sizes einsum with
    `optimize=True` re-derives its contraction path per call (the same trap
    `chain.py`'s own fold documents), and with `optimize=False` a
    three-operand contraction falls off BLAS entirely."""
    xp = bk.xp()
    B = E.shape[0]
    cl, d, cr = A_op.shape
    cl2 = A_conj.shape[0]
    crb = A_conj.shape[2]
    # X[b,j,p,k] = sum_i E[b,i,j] A_op[i,p,k]
    lhs = E.transpose(0, 2, 1).reshape(B * cl2, cl)
    X = lhs @ A_op.reshape(cl, d * cr)
    # E'[b,k,l] = sum_{j,p} X[b,j,p,k] conj(A)[j,p,l]
    X = X.reshape(B, cl2, d, cr).transpose(0, 3, 1, 2).reshape(B * cr, cl2 * d)
    out = X @ A_conj.reshape(cl2 * d, crb)
    return out.reshape(B, cr, crb)


# Deliberately NOT routed through `backend.jit`. The composites that
# module fuses all have shapes fixed by the bond dimension, which
# `set_pad_bonds` can freeze; this one's leading axis is the *environment
# batch*, which changes at every site and every trie node by construction.
# `jax.jit` retraces per input shape, so jitting here would trade one
# dispatch for thousands of traces -- the exact failure mode
# docs/pyitensor_gpu_port_plan.md Sec. 9 item 1 describes, with no padding
# available to remove it. The batch is what replaces the fusion: one GEMM
# here already carries what would otherwise be tens of thousands of calls.

# Cap on the elements of `_transfer`'s intermediate `X`, which is the
# largest array the sweep ever holds: `B * chi_r * chi_l * d`, i.e. `d *
# chi` times the environment batch itself. At n=30 the level-3 nodes carry
# ~22k environments, so an unchunked chi=80 step would allocate ~4.5 GB for
# one temporary. Splitting the batch costs nothing (the GEMMs stay large --
# a chunk is still tens of thousands of rows) and bounds the peak instead.
_CHUNK_ELEMS = 1 << 25   # 32M complex = 512 MB

# Byte budget for the live environments of one block of the sweep; see
# `_block_width`. Deliberately generous: a narrower block costs device
# dispatches, and the environments are the only large allocation here apart
# from the output tensor itself (which is n^4 complex, 1.6 GB at n=100, and
# is the caller's problem rather than this budget's).
_LIVE_BYTES = 2 << 30    # 2 GB


def _transfer_batched(E, A_op, A_conj, transfer):
    B = E.shape[0]
    per = A_op.shape[1] * A_op.shape[2] * A_conj.shape[0]
    step = max(1, _CHUNK_ELEMS // max(1, per))
    if B <= step:
        return transfer(E, A_op, A_conj)
    xp = bk.xp()
    return xp.concatenate([transfer(E[o:o + step], A_op, A_conj)
                           for o in range(0, B, step)], axis=0)


def _apply_op(A, mat):
    """The ket site tensor with a local operator folded in, `mat` indexed
    `[ket, bra]` exactly as `SiteSet.op()` builds it (`ITensor((s, s'),
    mat)`), i.e. the same call the fold makes."""
    xp = bk.xp()
    return xp.tensordot(A, mat, axes=([1], [0])).transpose(0, 2, 1)


class _MatrixRegistry:
    """Interns the distinct local matrices a shape enumeration produces, so
    trie nodes can be keyed on a small integer instead of an array."""

    def __init__(self):
        self._by_key = {}
        self.mats = []

    def intern(self, mat):
        """Returns an integer id, or None for an identically-zero matrix --
        a shape that produces one contributes nothing and is dropped before
        it ever reaches the sweep. On spinless sites this is what removes
        every shape placing two `Cdag`s (or two `C`s) on the same site,
        which is most of the repeated-index enumeration."""
        arr = np.ascontiguousarray(np.asarray(mat, dtype=complex))
        if not np.any(arr):
            return None
        key = arr.tobytes()
        got = self._by_key.get(key)
        if got is None:
            got = len(self.mats)
            self._by_key[key] = got
            self.mats.append(arr)
        return got


def _surjections(m):
    """Every assignment of the four factor positions onto exactly `m`
    ranks. Enumerating shapes this way rather than enumerating tuples is
    what makes the trie site-independent: `m! * S(4,m)` shapes (1, 14, 36,
    24 for m = 1..4) stand in for all n^4 tuples."""
    for assign in itertools.product(range(m), repeat=4):
        if len(set(assign)) == m:
            yield assign


def _shape_program(assign, names, site_matrix, F):
    """The (matrix id, parity) sequence one shape contracts, or None if any
    rank's composed matrix vanishes.

    `parity` is the number of factors placed up to and including this rank,
    mod 2 -- both the `F` folded into this rank's own matrix and the `F`
    applied to every gap site after it, which is the single rule the module
    docstring states."""
    m = max(assign) + 1
    steps = []
    placed = 0
    for r in range(m):
        here = [p for p in range(4) if assign[p] == r]
        mat = site_matrix(names[here[0]])
        for p in here[1:]:
            mat = site_matrix(names[p]) @ mat
        placed += len(here)
        parity = placed % 2
        if parity:
            mat = F @ mat
        steps.append((mat, parity))
    return steps


@lru_cache(maxsize=None)
def _combinations(n, m):
    """All strictly increasing m-tuples of 1-based sites, plus the mixed-
    radix code that indexes them, as arrays.

    The sweep builds the same code incrementally as the trie descends
    (`code * (n+1) + site`), never materializing the site tuples, so sorting
    a leaf node's codes is all it takes to line its values up with these
    combinations -- which `itertools.combinations` already emits in
    increasing-code order."""
    combos = np.array(list(itertools.combinations(range(1, n + 1), m)),
                      dtype=np.int64).reshape(-1, m)
    code = np.zeros(len(combos), dtype=np.int64)
    for r in range(m):
        code = code * (n + 1) + combos[:, r]
    return combos, code


def _block_width(trie, internal, n, arrays):
    """How many first-sites one block of the sweep may seed, from a byte
    budget on the live environments.

    A block seeding `k` first-sites keeps, for each internal trie node at
    level `r`, at most `k * C(n-1, r-1)` environments of `chi x chi`. Solving
    that for `k` is the whole calculation; the `chi` used is the largest bond
    in the MPS, so the estimate is an upper bound rather than a typical
    case."""
    chi = max(a.shape[0] for a in arrays)
    chi = max(chi, max(a.shape[2] for a in arrays))
    per_k = 0
    for key in internal:
        r = len(key)
        if r == 0:
            continue     # `internal` carries the trie root, which holds no
                         # environments of its own
        per_k += comb(n - 1, r - 1)
    per_k *= chi * chi * 16
    if per_k <= 0:
        return n
    return max(1, min(n, int(_LIVE_BYTES // per_k)))


def four_correlation_tensor_batched(chain, wf, cdag_ops, c_ops):
    """<Cdag_i C_j Cdag_k C_l> for every mode quadruple, batched.

    `cdag_ops`/`c_ops` are one `(operator_name, site_1based)` pair per mode,
    the same signature `Chain.four_correlation_tensor_fold` takes, so this
    covers a spinless chain (one mode per site) and a native spinful one
    (two modes per site) alike.

    Raises `ValueError` when the chain mixes site types, which would make a
    trie node's matrix depend on which site it landed on and so break the
    site-independence the whole enumeration rests on. Callers fall back to
    the fold in that case; no chain dmrgpy builds today hits it."""
    xp = bk.xp()
    nm = len(c_ops)
    if len(cdag_ops) != nm:
        raise ValueError("four_correlation_tensor_batched: cdag_ops and "
                          "c_ops must describe the same modes")
    n = chain.sites.length()
    types = {chain.sites.site_type(s) for s in range(1, n + 1)}
    if len(types) != 1:
        raise ValueError("four_correlation_tensor_batched: needs a uniform "
                          "site type across the chain")
    stype = types.pop()

    psi = wf.copy()
    psi.position(1)
    arrays = _arrays_lpr(psi)
    conj_arrays = [a.conj() for a in arrays]

    site_matrix = stype.matrix
    F = site_matrix("F")
    reg = _MatrixRegistry()
    # interned up front: the gap sites of an odd-parity trie node carry a
    # bare F, which no shape's own composed matrix is guaranteed to equal
    f_id = reg.intern(F)

    # ---- shapes -> trie -------------------------------------------------
    cdag_names = sorted({nm_ for nm_, _s in cdag_ops})
    c_names = sorted({nm_ for nm_, _s in c_ops})
    trie = {}          # node key -> {step: child key}
    leaf_shapes = {}   # node key -> [(assign, names, sign), ...]
    for n0 in cdag_names:
        for n1 in c_names:
            for n2 in cdag_names:
                for n3 in c_names:
                    names = (n0, n1, n2, n3)
                    for m in (1, 2, 3, 4):
                        for assign in _surjections(m):
                            steps = _shape_program(assign, names, site_matrix, F)
                            ids = [reg.intern(mat) for mat, _p in steps]
                            if any(i is None for i in ids):
                                continue  # a vanishing local product
                            key = ()
                            for r in range(m):
                                step = (ids[r], steps[r][1])
                                trie.setdefault(key, {})[step] = key + (step,)
                                key = key + (step,)
                            # sign of reordering (s_assign[0],..,s_assign[3])
                            # into site order: ranks are increasing sites, so
                            # this depends on the shape alone
                            sign = 1
                            for x in range(4):
                                for y in range(x + 1, 4):
                                    if assign[x] > assign[y]:
                                        sign = -sign
                            leaf_shapes.setdefault(key, []).append((assign, names, sign))
    mats = [xp.asarray(mm) for mm in reg.mats]
    internal = set(trie.keys())

    # ---- the sweep ------------------------------------------------------
    # Blocked on the first occupied site. Tuples whose smallest site falls in
    # a given range are swept together, and that is exact rather than
    # approximate: the trie merges two tuples only when their (matrix,
    # parity, site) prefixes agree, so environments with different first
    # sites never meet and splitting on it costs no arithmetic at all.
    #
    # What it buys is memory, and the block *width* is the knob between the
    # two things that can go wrong. One block (width n) holds 6*C(n,3)
    # level-3 environments -- 6.0 GB at n=100, chi=20, which no longer fits;
    # one block per site holds 6*C(n-a,2), 0.19 GB, but issues n times as
    # many array operations, and on a device each of those pays the ~0.35 ms
    # dispatch floor (n=100 would be ~150k dispatches, ~50 s of pure
    # dispatch -- the very cost this kernel exists to avoid). So the width is
    # chosen from a byte budget: as wide as `_LIVE_BYTES` allows, which is
    # one block for everything up to about n=40 and degrades gracefully
    # above it.
    def transfer(E, A_op, A_conj):
        return _transfer_batched(E, A_op, A_conj, _transfer)

    op_cache = {}

    def A_for(site, mat_id):
        """Site tensor with a local operator folded in, memoized: each
        (site, matrix) pair is now revisited once per block, so building it
        per visit would be n times the work it was before blocking."""
        key = (site, mat_id)
        got = op_cache.get(key)
        if got is None:
            A = arrays[site - 1]
            got = A if mat_id is None else _apply_op(A, mats[mat_id])
            op_cache[key] = got
        return got

    # operator-free left environments, one per site boundary: E0[s] contracts
    # sites 1..s with no operator anywhere. Built once and indexed by block,
    # replacing the per-tuple `psi.position(mn)` the fold uses.
    E0 = [xp.eye(1, dtype=complex)[None, :, :]]
    for s in range(1, n + 1):
        E0.append(transfer(E0[-1], arrays[s - 1], conj_arrays[s - 1]))

    leaves = {}    # node key -> [value arrays, code arrays]

    def record(ckey, Ec, code):
        shapes = leaf_shapes.get(ckey)
        if shapes is not None:
            slot = leaves.setdefault(ckey, [[], []])
            slot[0].append(xp.trace(Ec, axis1=1, axis2=2))
            slot[1].append(code)
        return ckey in internal

    width = _block_width(trie, internal, n, arrays)
    for a0 in range(1, n + 1, width):
        a1 = min(a0 + width, n + 1)     # first sites seeded by this block
        live = {}
        for s in range(a0, n + 1):
            Ac = conj_arrays[s - 1]
            born = []
            if s < a1:
                for step, ckey in trie[()].items():
                    Ec = transfer(E0[s - 1], A_for(s, step[0]), Ac)
                    born.append((ckey, Ec, np.array([s], dtype=np.int64)))
            for key, kids in trie.items():
                if key == ():
                    continue
                got = live.get(key)
                if got is None:
                    continue
                E, code = got
                for step, ckey in kids.items():
                    born.append((ckey, transfer(E, A_for(s, step[0]), Ac),
                                 code * (n + 1) + s))
            # gap propagation of everything already alive, before the newly
            # born join it -- a tuple's rank r+1 site is strictly right of its
            # rank r site, so a child created at s must not also be
            # transferred through s
            for key in list(live):
                E, code = live[key]
                gap = None if key[-1][1] == 0 else f_id
                live[key] = [transfer(E, A_for(s, gap), Ac), code]
            for ckey, Ec, code in born:
                if record(ckey, Ec, code):
                    got = live.get(ckey)
                    if got is None:
                        live[ckey] = [Ec, code]
                    else:
                        got[0] = xp.concatenate([got[0], Ec], axis=0)
                        got[1] = np.concatenate([got[1], code])

    # ---- scatter --------------------------------------------------------
    out = np.zeros((nm, nm, nm, nm), dtype=complex)
    lookup = {}
    for name in set(x for x, _ in cdag_ops) | set(x for x, _ in c_ops):
        lookup[name] = np.full(n + 1, -1, dtype=np.int64)
    for mode, (name, site) in enumerate(cdag_ops):
        lookup[name][site] = mode
    for mode, (name, site) in enumerate(c_ops):
        lookup[name][site] = mode

    for ckey, (vals, codes) in leaves.items():
        val = bk.to_host(xp.concatenate(vals) if len(vals) > 1 else vals[0])
        code = np.concatenate(codes) if len(codes) > 1 else codes[0]
        m = len(ckey)
        combos, _ = _combinations(n, m)
        # A leaf node is reached by exactly one increasing site tuple per
        # combination, so sorting its codes lines it up with `combos`, which
        # itertools.combinations already emits in that order. This replaces a
        # dense (n+1)**m lookup table -- 1.66 GB per m=4 node at n=100.
        if len(code) != len(combos):
            raise RuntimeError("four_correlation_tensor_batched: leaf node "
                                "%r produced %d values for %d site "
                                "combinations" % (ckey, len(code), len(combos)))
        got = val[np.argsort(code, kind="stable")]
        for assign, names, sign in leaf_shapes[ckey]:
            idx = [lookup[names[p]][combos[:, assign[p]]] for p in range(4)]
            ok = (idx[0] >= 0) & (idx[1] >= 0) & (idx[2] >= 0) & (idx[3] >= 0)
            if not ok.all():
                idx = [x[ok] for x in idx]
                out[idx[0], idx[1], idx[2], idx[3]] = sign * got[ok]
            else:
                out[idx[0], idx[1], idx[2], idx[3]] = sign * got
    return out
