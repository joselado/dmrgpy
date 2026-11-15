"""Free functions over MPS/MPO chains: sum, applyMPO, nmultMPO, inner
(=innerC -- see the note below), traceC, randomMPS. Mirrors the subset of
ITensor v3's own free-function API mpscpp3/chain_session.h calls directly.

Every one of these that grows bond dimension (sum, applyMPO, nmultMPO)
builds the exact result first and only then compresses it, with a single
truncating SVD sweep (_Chain.position) down to the requested
cutoff/maxdim, run in the gauge where discarding the smallest singular
values is the 2-norm-optimal truncation -- the standard, correct way to
compress a tensor-train sum/product, though not necessarily bit-for-bit
the same intermediate bond dimensions ITensor's own variational methods
would produce. (Truncating *while* building, which applyMPO/nmultMPO used
to do, is not the same thing and is not bounded by cutoff/maxdim at all;
see _apply_chain's docstring.) That's an intentional simplification (see autompo.py's
module docstring for the same reasoning applied to MPO construction):
dmrgpy only ever observes final numerical results bounded by
Cutoff/MaxDim, never internal bond dimensions. The one opt-in exception
is applyMPO(method="densitymatrix"), which truncates while it builds but
keeps the same basis at every bond (see _apply_densitymatrix).

ITensor v3 itself has both inner() (throws on complex operands) and
innerC() (always works) because IQTensor/ITensor's real-vs-complex type
split survives into v3 (see chain_session.h's long comment on this).
Phase 1 dropped that whole distinction on purpose -- every ITensor here is
unconditionally complex128 -- so there is exactly one implementation below,
exposed under both names purely so a later phase's transcription of
`innerC(...)` call sites needs no renaming.
"""

import numpy as np

from . import backend as bk

from .index import Index
from .mpscontainer import MPO, MPS, _link_at
from .svd import eigh_truncate, qr_split, svd
from .tensor import ITensor, commonIndex, contract_many
from .tensor import dag as _t_dag
from .tensor import prime as _t_prime


def _fresh_link_copy(chain):
    """A new list of tensors representing the same chain, but with every
    Link index replaced by a freshly minted one (physical indices
    untouched, and consistent across each bond so the chain still hangs
    together). Needed for the *bra* side of every inner()/traceC-adjacent
    contraction below: dag() conjugates values but leaves Index identity
    alone, so `dag(X.A(i)) * X.A(i)` for literally the same chain object X
    (e.g. inner(psi,psi), the single most common call pattern in
    chain_session.h -- every normalization and every same_mps() check is a
    self-overlap) would otherwise auto-contract the *Link* legs too, not
    just the physical one, silently collapsing each site to its own
    Frobenius norm instead of building the correct environment. Relabeling
    the bra's links first means only the physical leg can ever match the
    ket's (or an operator's) indices, exactly as intended."""
    n = chain.length()
    old_links = {k: _link_at(chain, k, k + 1) for k in range(1, n)}
    new_links = {k: old_links[k].sim() for k in range(1, n)}
    tensors = []
    for i in range(1, n + 1):
        T = chain.A(i)
        old_l, old_r = old_links.get(i - 1), old_links.get(i)
        new_inds = []
        for ind in T.inds:
            if old_l is not None and ind == old_l:
                new_inds.append(new_links[i - 1])
            elif old_r is not None and ind == old_r:
                new_inds.append(new_links[i])
            else:
                new_inds.append(ind)
        tensors.append(ITensor(tuple(new_inds), T.array))
    return tensors


def inner(*args):
    """inner(A,B) = <A|B> for two MPS, or inner(bra,M,ket) = <bra|M|ket>
    for an MPS/MPO/MPS triple. Each site's local contribution and the
    running environment are all contracted together via contract_many()
    (see tensor.py), not a naive left-to-right `*` chain -- correct
    regardless of gauge/canonical form either way, since it's just one big
    contraction performed incrementally, but the 3-arg case in particular
    used to build `piece = bra_i*M.A(i)*ket.A(i)` *before* ever touching
    the accumulated environment, which is exactly the same
    intermediate-tensor blowup dmrg.py's environment extension had (see
    contract_many()'s docstring for the measured numbers) -- confirmed
    directly here too: inner(psi,H,psi) on a 14-site, maxdim=60 chain,
    14.09s before this fix, 0.007s after, same result."""
    if len(args) == 2:
        A, B = args
        n = A.length()
        bra_tensors = _fresh_link_copy(A)
        env = None
        for i in range(1, n + 1):
            piece = _dag_local(bra_tensors[i - 1])
            pieces = [p for p in (env, piece, B.A(i)) if p is not None]
            env = contract_many(pieces)
        return env.scalar()
    if len(args) == 3:
        bra, M, ket = args
        n = bra.length()
        bra_tensors = _fresh_link_copy(bra)
        env = None
        for i in range(1, n + 1):
            # M.A(i)'s "in" (unprimed) leg must contract against ket's
            # physical leg, and its "out" (primed) leg against the bra's --
            # so the bra's physical leg needs priming first, exactly like
            # dmrg.py's own <psi|H|psi> environments (_relabel_bra_local)
            # do. Without this, both legs accidentally match the *ket*
            # (same plev), leaving M's output leg and the bra's own
            # physical leg dangling uncontracted instead.
            piece = _dag_local(_t_prime(bra_tensors[i - 1], "Site"))
            pieces = [p for p in (env, piece, M.A(i), ket.A(i)) if p is not None]
            env = contract_many(pieces)
        return env.scalar()
    raise ValueError("inner: expected 2 or 3 arguments, got {}".format(len(args)))


innerC = inner


def _mps_arrays_lpr_device(psi):
    """Every site tensor of `psi` as a (chi_left, d, chi_right) array in
    the *active* array backend, length-1 axes inserted at the open ends.

    chain.py has a `_mps_arrays_lpr` that does the same reordering and then
    calls `np.ascontiguousarray`, which is correct there (the four-point
    kernel converts once, deliberately, and wants a host array) and exactly
    wrong here: this feeds a batch that must stay wherever the MPS already
    lives. Same explicit transpose-against-named-indices, since index order
    on an MPS tensor is not guaranteed."""
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
        out.append(arr)
    return out


def _batched_overlap_site(E, bra, ket):
    """One site of BatchedBras.overlaps(): E[b,aL,cL] -> E'[b,aR,cR], as
    the two batched GEMMs described in that class's docstring.

    bra is already conjugated and padded, shape (B, AL, d, AR); ket is
    (cL, d, cR); E is (B, AL, cL)."""
    _xp = bk.xp()
    B, AL, d, AR = bra.shape
    cL, _, cR = ket.shape
    # X[b, cL, d*AR] = sum_aL E[b, aL, cL] * bra[b, aL, d*AR]
    X = _xp.matmul(E.transpose(0, 2, 1), bra.reshape(B, AL, d * AR))
    # E'[b, AR, cR] = sum_{cL,d} X[b, cL, d, AR] * ket[cL, d, cR]
    X = X.reshape(B, cL * d, AR)
    return _xp.matmul(X.transpose(0, 2, 1), ket.reshape(cL * d, cR))


_batched_overlap_site = bk.jit(_batched_overlap_site)


class BatchedBras:
    """A fixed set of bra MPS, prepared once so that <bra_b|ket> for every
    b comes from a single sweep instead of one sweep per bra.

    Why this exists. tdz.py's complex-time correlator (submode "TDZ",
    arXiv:2311.10909) needs n_max+1 overlaps of the *same* evolving state
    against a *fixed* family of bras -- phi^(n)(t) = <H^n B GS|psi(t)> --
    at every one of hundreds of time steps. Written as a loop over
    `inner()` that is (n_max+1) independent left-to-right sweeps, i.e.
    (n_max+1) * N environment updates of a few array operations each, all
    of them chi x chi and none of them big enough to occupy a GPU. Batched,
    it is N steps of two GEMMs carrying all n_max+1 environments at once:
    the same arithmetic, a fifth of the dispatches, which on a device is
    the quantity that decides the cost (see docs/pyitensor_gpu_port_plan.md
    Sec. 9's four-point entry for the same lesson learned the expensive
    way).

    The bras never change, so the transpose-to-(left,site,right), the
    conjugation and the zero-padding to a common per-bond width all happen
    here, once, rather than per step. Zero padding is exact rather than an
    approximation: a padded column of the bra contracts against the ket to
    zero and contributes nothing to any overlap, which is the same argument
    backend.set_pad_bonds relies on.

    Bond widths are equalized *per bond*, not to one global maximum: near
    the chain ends the true dimensions are 2, 4, 8, ... and padding those up
    to the middle's width would do a great deal of arithmetic on blocks
    known to be zero.

    The prepared copy roughly doubles what the bras already occupy (B x N
    padded site tensors: ~276 MB for B=5, N=30 at chi=240, complex128).
    That is the cost of preparing once instead of per step, and it is
    bounded by the caller's own maxm, but it is worth knowing before
    batching a large family of bras rather than a handful.
    """

    __slots__ = ("_bras", "_n", "_nsites", "_batch")

    def __init__(self, chains):
        if not chains:
            raise ValueError("BatchedBras: needs at least one bra")
        _xp = bk.xp()
        self._batch = len(chains)
        self._nsites = chains[0].length()
        for c in chains:
            if c.length() != self._nsites:
                raise ValueError("BatchedBras: every bra must have the same "
                                 "number of sites")
        per_bra = [_mps_arrays_lpr_device(c) for c in chains]
        stacked = []
        for i in range(self._nsites):
            mats = [arrs[i] for arrs in per_bra]
            AL = max(m.shape[0] for m in mats)
            d = mats[0].shape[1]
            AR = max(m.shape[2] for m in mats)
            padded = []
            for m in mats:
                if m.shape != (AL, d, AR):
                    buf = bk.zeros((AL, d, AR))
                    buf = bk.setblock(buf, (slice(0, m.shape[0]), slice(None),
                                            slice(0, m.shape[2])), m)
                    m = buf
                padded.append(m.conj())
            stacked.append(_xp.stack(padded))
        self._bras = stacked
        self._n = self._nsites

    def __len__(self):
        return self._batch

    def overlaps(self, ket):
        """[<bra_b|ket> for b] as a NumPy complex array of length B, from
        one batched sweep. One host transfer at the end, of B numbers."""
        _xp = bk.xp()
        kets = _mps_arrays_lpr_device(ket)
        if len(kets) != self._nsites:
            raise ValueError("BatchedBras.overlaps: ket has %d sites, bras "
                             "have %d" % (len(kets), self._nsites))
        E = _xp.ones((self._batch, 1, 1), dtype=complex)
        for i in range(self._nsites):
            E = _batched_overlap_site(E, self._bras[i], kets[i])
        return np.asarray(bk.to_host(E.reshape(self._batch)))

    def sandwiches(self, mpo, ket):
        """[<bra_b|mpo|ket> for b] from one batched sweep, `mpo` given as
        mpo_arrays_device(W) so that a caller sandwiching many kets with
        one operator converts it once. The same numbers as a loop of
        inner(bra_b, W, ket), with every bra carried through each GEMM at
        once; on small chains the loop is dominated by Python overhead
        (2729 inner calls were 65% of a 4-site Lindblad spectrum)."""
        _xp = bk.xp()
        kets = _mps_arrays_lpr_device(ket)
        if len(kets) != self._nsites or len(mpo) != self._nsites:
            raise ValueError("BatchedBras.sandwiches: the ket, the MPO and "
                             "the bras need the same number of sites")
        E = _xp.ones((self._batch, 1, 1, 1), dtype=complex)
        for i in range(self._nsites):
            E = _batched_sandwich_site(E, self._bras[i], mpo[i], kets[i])
        return np.asarray(bk.to_host(E.reshape(self._batch)))


def mpo_arrays_device(W):
    """Every site tensor of the MPO W as a (w_left, s_out, s_in, w_right)
    array in the active backend, length-1 axes at the open ends; s_out is
    the primed (output) site index, contracted with the bra."""
    n = W.length()
    out = []
    for i in range(1, n + 1):
        T = W.A(i)
        sites = [ind for ind in T.inds if ind.hastags("Site")]
        s_in = next(ind for ind in sites if ind.plev == 0)
        s_out = next(ind for ind in sites if ind.plev == 1)
        left = commonIndex(W.A(i - 1), T) if i > 1 else None
        right = commonIndex(T, W.A(i + 1)) if i < n else None
        order = [ind for ind in (left, s_out, s_in, right) if ind is not None]
        arr = T.transpose_to(order)
        if left is None:
            arr = arr.reshape((1,) + arr.shape)
        if right is None:
            arr = arr.reshape(arr.shape + (1,))
        out.append(bk.asarray(arr))
    return out


def _batched_sandwich_site(E, bra, W, ket):
    """One site of BatchedBras.sandwiches(): E[b,a,w,c] -> E'[b,A,v,C],
    with bra (B, a, e, A) conjugated and padded, W (w, e, d, v) and ket
    (c, d, C): three GEMMs, ket first, then the MPO, then the bras."""
    _xp = bk.xp()
    B, a, w, c = E.shape
    _, d, C = ket.shape
    _, e, _, v = W.shape
    A = bra.shape[3]
    X = _xp.matmul(E.reshape(B * a * w, c), ket.reshape(c, d * C))
    X = X.reshape(B, a, w, d, C).transpose(0, 1, 4, 2, 3)
    Y = _xp.matmul(X.reshape(B * a * C, w * d),
                   W.transpose(0, 2, 1, 3).reshape(w * d, e * v))
    Y = Y.reshape(B, a, C, e, v).transpose(0, 2, 4, 1, 3)
    Z = _xp.matmul(Y.reshape(B, C * v, a * e), bra.reshape(B, a * e, A))
    return Z.reshape(B, C, v, A).transpose(0, 3, 2, 1)


_batched_sandwich_site = bk.jit(_batched_sandwich_site)


def _dag_local(T):
    from .tensor import dag
    return dag(T)


def traceC(mpo):
    """Tr[A] for an MPO -- mirrors v3's traceC(), which mpscpp3/
    chain_session.h's trace_operator() uses directly instead of building an
    explicit identity MPO the way the old file-based backend had to."""
    n = mpo.length()
    env = None
    for i in range(1, n + 1):
        T = mpo.A(i)
        site_inds = [ind for ind in T.inds if ind.hastags("Site")]
        in_ind = next(ind for ind in site_inds if ind.plev == 0)
        out_ind = next(ind for ind in site_inds if ind.plev == 1)
        ax_in = T.inds.index(in_ind)
        ax_out = T.inds.index(out_ind)
        remaining = [ind for ind in T.inds if ind != in_ind and ind != out_ind]
        arr = bk.xp().trace(T.array, axis1=ax_in, axis2=ax_out)
        piece = ITensor(tuple(remaining), arr)
        env = piece if env is None else env * piece
    return env.scalar()


def mpo_dagger(A):
    """A^dagger for an MPO: every tensor conjugated, with its two Site
    legs swapped."""
    from .tensor import dag, swapPrime
    out = A.copy()
    for j in range(1, out.length() + 1):
        out.set_A(j, dag(swapPrime(out.A(j), 0, 1, "Site")))
    return out


def _mpo_frobenius_norm(M):
    """||M||_F = sqrt(Tr M^dagger M), read off the orthogonality center
    after one lossless right-to-left sweep. The sweep is forced from the
    far end, since a chain fresh out of sum_many() claims center=n while
    being in no canonical form at all; from there every site it passes is
    right-orthogonalized whatever its gauge (_shift_right's argument in
    mpscontainer.py)."""
    out = M.copy()
    n = out.length()
    out.center = n
    out.position(1)
    return float(bk.to_host(bk.xp().linalg.norm(out.A(1).array)))


# Relative size of H - H^dagger below which an MPO counts as Hermitian.
# The norm of the difference is computed without cancellation (see
# mpo_is_hermitian), so its rounding floor is ~1e-15 times ||H||_F, and
# this sits three orders above it -- the level of
# multioperator.clean_threshold.
_HERMITIAN_RTOL = 1e-12


def mpo_is_hermitian(H, rtol=_HERMITIAN_RTOL):
    """Whether ||H - H^dagger||_F <= rtol * ||H||_F for an MPO H.

    The difference is built as an exact direct sum and canonicalized, so
    its norm comes from the center tensor rather than from 2 Tr(H^dag H) -
    2 Re Tr(H H), a difference of two numbers of order ||H||_F^2 that
    would only resolve the ratio down to ~1e-8. The cost is one SVD sweep
    over an MPO of twice H's bond dimension, far below one TDVP step.
    An exact zero operator counts as Hermitian."""
    scale = _mpo_frobenius_norm(H)
    if scale == 0.0:
        return True
    D = sum_many([H, mpo_dagger(H) * (-1.0)])
    return _mpo_frobenius_norm(D) <= rtol * scale


def _bond_dims(chain):
    """{k: dim of the Link between site k and k+1} for k=1..N-1."""
    n = chain.length()
    return {k: _link_at(chain, k, k + 1).dim for k in range(1, n)}


def sum_many(chains, cutoff=0.0, maxdim=None):
    """N-way direct sum of MPS (or MPO) chains: the same exact block-
    diagonal concatenation sum() does for two operands, generalized to K
    via one K-way placement per site instead of K-1 pairwise merges each
    paying their own truncating SVD sweep.

    This was written for mpobuilder.py's to_mpo(), which no longer uses
    it: building a T-term Hamiltonian by folding in one term at a time via
    T-1 calls to sum() forces T-1 full left-to-right SVD sweeps, even
    though that per-step compression demonstrably doesn't reduce bond
    dimension at all until a final bidirectional pass runs -- confirmed
    directly, a 14-site nearest-neighbor Heisenberg chain (39 terms)
    landed at bond dimension exactly 39, i.e. zero compression, before
    that final pass brought it down to 5. Concatenating all T operands at
    once is exact array placement (no linear algebra), so only the
    caller's own final compression sweep(s) need to touch SVD at all,
    turning that builder's O(T) SVD sweeps into O(1).

    That fix capped the number of sweeps but not their *size*: the
    concatenated MPO still had bond dimension T going into them, which is
    what made MPO construction O(L^4) until mpobuilder.py stopped
    concatenating terms at all (it assembles the finite-state machine
    directly now -- see its module docstring). sum_many() remains in use
    for genuine K-way sums elsewhere, e.g. sector.py's in-sector start
    state and charge-penalty operator."""
    chains = list(chains)
    if not chains:
        raise ValueError("sum_many: no chains given")
    if len(chains) == 1:
        return chains[0].copy()
    n = chains[0].length()
    for c in chains[1:]:
        if c.length() != n:
            raise ValueError("sum_many: mismatched chain length {} vs {}".format(n, c.length()))

    new_links = {}
    for k in range(1, n):
        total = 0
        for c in chains:
            total += _link_at(c, k, k + 1).dim
        new_links[k] = Index(total, tags="Link,l={}".format(k))

    tensors = []
    for i in range(1, n + 1):
        Ts = [c.A(i) for c in chains]
        phys = tuple(ind for ind in Ts[0].inds if ind.hastags("Site"))
        phys_set = set(phys)
        for T in Ts[1:]:
            if set(ind for ind in T.inds if ind.hastags("Site")) != phys_set:
                raise ValueError("sum_many: different index structure at site {}".format(i))

        links_l = [_link_at(c, i, i - 1) for c in chains]
        links_r = [_link_at(c, i, i + 1) for c in chains]

        new_left = new_links.get(i - 1)
        new_right = new_links.get(i)
        phys_shape = tuple(p.dim for p in phys)
        shape = (((new_left.dim,) if new_left else ())
                 + phys_shape
                 + ((new_right.dim,) if new_right else ()))
        combined = bk.zeros(shape)

        left_off = 0
        right_off = 0
        for T, la, ra in zip(Ts, links_l, links_r):
            order = ([la] if la else []) + list(phys) + ([ra] if ra else [])
            arr = T.transpose_to(order)
            left_dim = la.dim if la else 0
            right_dim = ra.dim if ra else 0
            idx = []
            if new_left:
                idx.append(slice(left_off, left_off + left_dim))
            idx += [slice(None)] * len(phys)
            if new_right:
                idx.append(slice(right_off, right_off + right_dim))
            # setblock, not `combined[...] = arr`: JAX arrays are immutable,
            # so this returns a new array rather than mutating in place.
            combined = bk.setblock(combined, tuple(idx), arr)
            left_off += left_dim
            right_off += right_dim

        inds = (([new_left] if new_left else []) + list(phys)
                + ([new_right] if new_right else []))
        tensors.append(ITensor(tuple(inds), combined))

    cls = type(chains[0])
    result = cls(tensors)
    result.center = n
    if cutoff or maxdim is not None:
        # The block-diagonal concatenation above is in no canonical form at
        # all, so a *truncating* left-to-right sweep straight off it would
        # weigh each cut's singular values by whatever norm the untouched
        # right part carries instead of by the state's true Schmidt values,
        # and discard the wrong components. Right-canonicalize losslessly
        # first (cutoff=0/maxdim=None, gauge only), then truncate on the way
        # back -- the standard MPS compression recipe, and the same fix
        # _apply_chain() needs for the same reason; see its comment for the
        # measured symptom. Skipped entirely when nothing can be discarded
        # (cutoff=0 and no maxdim), where the single sweep is already an
        # exact canonicalization.
        result.position(1)
    result.position(n, cutoff=cutoff, maxdim=maxdim)
    return result


def sum(A, B, cutoff=0.0, maxdim=None):
    """Direct sum of two MPS (or two MPO): the exact, standard MPS/MPO
    addition construction (concatenation at the boundaries, block-diagonal
    in the link spaces at interior sites), followed by one truncating
    left-to-right SVD sweep down to cutoff/maxdim. A thin 2-operand
    wrapper around sum_many()."""
    return sum_many([A, B], cutoff=cutoff, maxdim=maxdim)


def _apply_chain(K, X, out_cls, cutoff=0.0, maxdim=None):
    """Shared implementation of applyMPO (X=MPS, out_cls=MPS) and nmultMPO
    (X=MPO, out_cls=MPO): a single left-to-right sweep that contracts K
    against X one site at a time (whatever physical legs match,
    auto-contract per ITensor.__mul__) and splits each cut before moving
    on, carrying the remainder -- `leftover`, with legs (new bond, K's own
    right link, X's own right link) -- forward into the next site's
    contraction. That sweep is *exact* -- it truncates nothing, and uses a
    QR rather than an SVD since no singular values are wanted -- and the
    caller's cutoff/maxdim are then enforced by one truncating
    right-to-left sweep over the finished chain.

    Truncating inside the sweep instead -- the "zip-up" of Stoudenmire &
    White, New J. Phys. 12, 055026 (2010), §3.2, which is what this
    function used to do -- is *not* equivalent, and the difference is not
    a tolerance. At cut i the left side of the SVD is orthonormal (it is
    the previous U), but the right side is the not-yet-contracted product
    of K's and X's remaining tensors, and that product is NOT
    right-orthogonal even when K and X separately are: MPO/MPS
    right-canonicity gives sum_{s,s',aR} K K* = delta_{aL aL'}, whereas
    the product needs the finer sum_{s',aR} K K* = delta_{aL aL'}
    delta_{s s2}, which nothing supplies. So the singular values of
    `piece` are not the Schmidt values of the exact product, and dropping
    the smallest of them is not the 2-norm-optimal truncation. Measured
    (2026-09 audit, finding 13): applying an 8-site Heisenberg H to its
    own ground state, whose exact H|psi> has Schmidt rank [2,4,8,16,8,4,2]
    -- so maxdim=16 discards literally nothing, the 17th singular value
    being 0.0 at every bond -- the zip-up still lost 9.4e-05 in 2-norm,
    against 1.6e-14 at maxdim=32. Downstream that is a
    gs_energy_fluctuation() of 5.1e-05 on a state that is an eigenstate,
    i.e. a user watching that number to tune maxm is reading the
    application error, not the state's.

    The exact sweep costs more (its cut dimensions run up to
    dim(K-link)*dim(X-link) rather than maxdim), and that is the price of
    the truncation meaning what it says. Measured end to end on a KPM
    dynamical correlator (12 sites, maxm=kpmmaxm=30: 1.30s -> 1.82s; 16
    sites, maxm=kpmmaxm=60: 9.9s -> 13.5s), i.e. ~1.4x, of which the QR
    (rather than a second SVD) in the exact sweep buys back a factor of
    ~2. The final sweep is skipped entirely when there is nothing to
    enforce (cutoff=0 and no maxdim), where the exact sweep already *is*
    the whole answer -- note that in that case the result keeps QR-sized
    bonds rather than rank-revealed ones, an identical operator stored
    slightly more loosely; every call site inside this package passes a
    maxdim, so that path is library-use only.

    This differs from the textbook-simplest approach (contract K.A(i)*X.A(i)
    at *every* site first, mechanically fusing each site's K-link and X-link
    into one combined dim(K-link)*dim(X-link) Link index, only *then*
    running one truncating sweep over the whole chain) only in when the
    left side of that fused index gets collapsed back down to
    cutoff/maxdim -- both reach the same fixed point (a left-to-right
    canonical sweep absorbs the compressed remainder of site i into site
    i+1 before site i+1 is ever finalized either way), but the
    every-site-first approach pays to build a (dim(K-link)*dim(X-link),
    phys, dim(K-link)*dim(X-link)) tensor at *every* site via a full
    tensordot, including sites whose left side is about to be collapsed
    from a full dim(K-link)*dim(X-link) back down to at most maxdim by the
    very next SVD -- confirmed directly, that discarded left-side
    tensordot work was ~15-20% of this function's total time on a
    representative KPM dynamical-correlator profile (8-site chain,
    maxm=kpmmaxm=30). Interleaving construction with truncation instead
    means every site's own contraction only ever touches an already-
    truncated (<=maxdim) left side, never a fused-but-about-to-be-
    discarded one.
    """
    # Right-canonicalize both inputs before the left-to-right sweep. This
    # keeps the sweep's cut dimensions honest (an SVD at cutoff=0 drops
    # only genuinely-zero singular values when the right side carries no
    # spurious norm of its own) and it is the canonicalization the
    # zip-up algorithm assumes (Stoudenmire & White, New J. Phys. 12,
    # 055026 (2010), §3.2). It was originally added because, without it,
    # applying a Heisenberg (H - E0 - omega) MPO to an MPS on a 10-site
    # chain -- bond dimension 32, far *below* maxdim 60, so maxdim never
    # even bound -- made <Hb|Hb> and <b|H^2 b> disagree by 0.86%
    # (0.0018381 vs 0.0018222) where the compiled ITensor backends agree
    # to 1e-15, and left cvm.py's conjugate gradient stalling at a
    # residual of ~7e-2 on itensor_version="python" (against ~9e-6 on
    # itensor_version 2 and 3) with a visibly wrong spectrum. The comment
    # here used to go further and claim that this made the per-cut
    # truncation *optimal* ("sum (K†K) (x) (X†X) = I (x) I"); it does not
    # -- see this function's docstring for why the product of two
    # right-canonical chains is not itself right-canonical, and for the
    # measurement that finally pinned it. That is why the sweep below no
    # longer truncates at all.
    #
    # Done in place rather than on copies: position() with the default
    # cutoff=0/maxdim=None is lossless and only changes the gauge, so the
    # represented state/operator is unchanged, and mutating in place means a
    # reused operator (cvm.py applies the same Hshift hundreds of times) pays
    # for it once -- every later call finds center==1 and both while loops in
    # position() are no-ops.
    for chain in (K, X):
        if chain.center is None:
            chain.center = chain.length()
        if chain.center != 1:
            chain.position(1)

    n = K.length()
    tensors = []
    leftover = None  # (bond, K's right link, X's right link) from the previous site, or None at site 1
    for i in range(1, n + 1):
        Kt, Xt = K.A(i), X.A(i)
        if leftover is not None:
            piece = leftover * Kt
            piece = piece * Xt
        else:
            piece = Kt * Xt
        kR, xR = _link_at(K, i, i + 1), _link_at(X, i, i + 1)
        if kR is None and xR is None:  # last site: nothing left to split off
            tensors.append(piece)
            break
        right_links = set(l for l in (kR, xR) if l is not None)
        left_inds = [ind for ind in piece.inds if ind not in right_links]
        # Exact split: cutoff/maxdim are deliberately NOT applied here,
        # they are enforced by the right-to-left sweep below instead (see
        # this function's docstring). QR rather than an untruncated SVD
        # since nothing here looks at the singular values -- same split,
        # about half the FLOPs, see svd.py's qr_split().
        Q, C, _bond = qr_split(piece, left_inds, orthonormal="left")
        tensors.append(Q)
        leftover = C

    result = out_cls(tensors)
    result.center = n
    if maxdim is not None or cutoff > 0.0:
        # Sites 1..n-1 are the U's of the sweep above, i.e. left-orthogonal,
        # so at every cut of this right-to-left sweep the left environment
        # is an isometry and the right one is already right-orthogonal from
        # the part of the sweep that has run -- which is exactly the gauge
        # in which discarding the smallest singular values IS the optimal
        # 2-norm truncation. Same recipe (and same reason) as sum_many()'s
        # own position(1)-then-position(n) above.
        result.position(1, cutoff=cutoff, maxdim=maxdim)
    return result


def _relabel(T, mapping):
    """T with every index that is a key of `mapping` replaced by its
    value, same data. Index equality includes the prime level, so a map
    keyed on an unprimed leg leaves its primed partner alone."""
    return ITensor(tuple(mapping.get(ind, ind) for ind in T.inds), T.array)


def _apply_densitymatrix(K, x, cutoff=0.0, maxdim=None):
    """K (MPO) applied to x (MPS) by the density-matrix algorithm of
    Stoudenmire and White, New J. Phys. 12, 055026 (2010), the default of
    ITensorMPS.jl's apply(A, psi; alg="densitymatrix"), whose mpo.jl this
    follows. Same output convention as applyMPO: the physical legs of the
    result are K's primed out legs.

    First the exact left environments E_j of <x|K^dag K|x> over sites
    1..j. Then a sweep from the right end: at site j the carry R (the
    product on sites j..n, already projected on the basis kept to the
    right of j) and E_{j-1} give the reduced density matrix rho of the
    product on (s'_j, right bond). Its leading eigenvectors, kept by
    svd.eigh_truncate's usual cutoff/maxdim rule (the eigenvalues are the
    squared singular values, so the rule is the same), become the new
    right-orthogonal site tensor, and R moves one site left. Site 1 takes
    the remainder: the result has center 1 and the norm of the projected
    product, not normalized.

    These eigenvectors are the right Schmidt vectors of the same partially
    projected product that _apply_chain's truncating SVD sweep keeps at
    the same bond. So the two methods give the same state, up to roundoff
    and to the choice inside a degenerate multiplet at a cut; this is not
    the zip-up of _apply_chain's docstring. What changes is the cost: the
    product never exists at bond dimension chi_K*chi_x, and a site costs
    about O(d chi_x^3 chi_K^2) against O(d (chi_K chi_x)^3), a factor of
    about chi_K less. The price is memory, n-1 environments of
    chi_x^2 chi_K^2 entries, and precision: rho's eigenvalues are known
    to about eps times the largest, so singular values below about 1e-8 of
    the largest are not resolved. At cutoff=0 a product with such a tail
    comes out good to about 1e-8 relative in 2-norm, not to roundoff; any
    cutoff of 1e-15 or more drops that tail in both methods.

    The kept dimension at a bond is also capped by dim(x link)*dim(K link)
    on its left, the largest rank the product can have there. Inside that
    cap cutoff=0 still keeps the null directions of rho whose eigenvalue
    came out as positive roundoff, so near site 1, where the cap can exceed
    the Hilbert-space bound, the bonds can come out larger than the exact
    route's; the state is the same, and a cutoff of about 1e-14 or more
    removes them. K and x are not
    modified and need no canonical form, since the environments are exact
    in any gauge. Every step is an ITensor contraction or an eigh, so this
    runs on the JAX backend as it is; only rho's spectrum, O(chi), comes to
    the host. set_pad_bonds is not honoured: the bonds come out at the kept
    dimension, so under jit the shapes follow it."""
    n = K.length()
    if x.length() != n:
        raise ValueError("_apply_densitymatrix: MPO has %d sites, MPS has %d"
                         % (n, x.length()))
    Ks = [K.A(i) for i in range(1, n + 1)]
    xs = [x.A(i) for i in range(1, n + 1)]
    if n == 1:
        out = MPS([Ks[0] * xs[0]])
        out.center = 1
        return out

    # 0-based lists: kl[j], xl[j] are the links between sites j and j+1.
    kl = [commonIndex(Ks[j], Ks[j + 1]) for j in range(n - 1)]
    xl = [commonIndex(xs[j], xs[j + 1]) for j in range(n - 1)]
    if any(l is None for l in kl + xl):
        raise ValueError("_apply_densitymatrix: every interior bond of K and "
                         "x needs a link")
    if set(kl) & set(xl):
        # A link shared by K and x would be contracted by K*x as if it
        # were a physical leg. Chains built separately never share one.
        Ks = _fresh_link_copy(K)
        kl = [commonIndex(Ks[j], Ks[j + 1]) for j in range(n - 1)]
    s_in = [commonIndex(Ks[j], xs[j], "Site") for j in range(n)]
    if any(s is None for s in s_in):
        raise ValueError("_apply_densitymatrix: K and x share no site index "
                         "at some site")
    s_out = [next(ind for ind in Ks[j].inds
                  if ind.hastags("Site") and ind != s_in[j])
             for j in range(n)]

    # The bra layer of <x|K^dag K|x>: every link and every in-leg replaced
    # by a fresh copy. K's out leg s'_j is not in the map (it differs from
    # s_in[j] by its prime level), so the bra K still contracts with the
    # ket K through it.
    bra = {ind: ind.sim() for ind in kl + xl + s_in}

    # E[j]: sites 0..j, indices (xl[j], kl[j], bra kl[j], bra xl[j]).
    E = []
    env = None
    for j in range(n - 1):
        Kc = _t_dag(_relabel(Ks[j], bra))
        xc = _t_dag(_relabel(xs[j], bra))
        env = xs[j] if env is None else env * xs[j]
        env = ((env * Ks[j]) * Kc) * xc
        E.append(env)

    tensors = [None] * n
    R = xs[n - 1] * Ks[n - 1]  # (xl[n-2], s'_{n-1}, kl[n-2])
    b = None                   # output bond on the right of site j
    for j in range(n - 1, 0, -1):
        # R: (xl[j-1], kl[j-1], s'_j, b). Its bra copy Rc is R conjugated
        # with every index renamed, rather than a second contraction chain.
        rows = [s_out[j]] + ([b] if b is not None else [])
        cols = [ind.sim() for ind in rows]
        ren = {xl[j - 1]: bra[xl[j - 1]], kl[j - 1]: bra[kl[j - 1]]}
        ren.update(zip(rows, cols))
        Rc = _t_dag(_relabel(R, ren))
        rho = (E[j - 1] * R) * Rc  # (rows, cols)
        dims = tuple(ind.dim for ind in rows)
        D = int(np.prod(dims))
        mat = rho.transpose_to(rows + cols).reshape(D, D)
        mat = 0.5 * (mat + mat.conj().T)
        cap = xl[j - 1].dim * kl[j - 1].dim
        if maxdim:
            cap = min(cap, int(maxdim))
        U, keep = eigh_truncate(mat, cutoff, cap)
        # rho = Phi Phi^dag with Phi = R as a (rows, rest) matrix, so the
        # product is approximated by U U^dag Phi: the site tensor is U as
        # it is, on (new left bond, s'_j, b), and the carry is U^dag Phi.
        bnew = Index(keep, tags="Link,l={}".format(j))
        arr = U.reshape(dims + (keep,)).transpose((len(dims),) + tuple(range(len(dims))))
        B = ITensor((bnew,) + tuple(rows), arr)
        tensors[j] = B
        R = ((R * _t_dag(B)) * xs[j - 1]) * Ks[j - 1]  # (bnew, xl[j-2], kl[j-2], s'_{j-1})
        b = bnew
    tensors[0] = R
    out = MPS(tensors)
    out.center = 1
    return out


# Which algorithm applyMPO uses when the caller names none, read at call
# time: "exact" (_apply_chain) or "densitymatrix" (_apply_densitymatrix).
APPLY_METHOD = "exact"


def applyMPO(K, x, x0=None, cutoff=0.0, maxdim=None, method=None):
    """K (MPO) applied to x (MPS). The physical leg of the result is left
    *primed* (matching K's own "out" leg) -- callers noPrime() it back to a
    plain ket afterward, exactly like every apply_mpo() call site in
    mpscpp3/chain_session.h does.

    `method` picks the algorithm: "exact" is _apply_chain, the product
    formed exactly and then truncated by one SVD sweep; "densitymatrix"
    is _apply_densitymatrix, the same truncation at a cost lower by about
    the MPO bond dimension. None takes the module default APPLY_METHOD,
    which is "exact"; any other name raises.

    `x0` (an initial guess, used by real ITensor to seed an iterative
    variational "Fit" solve) is accepted for call-signature compatibility
    with chain_session.h's two applyMPO() overloads but is otherwise
    unused: both methods form the truncation of the exact product, in the
    basis that makes it optimal, so the result is bounded by cutoff/maxdim
    independently of x0 and ignoring it is correctness-preserving."""
    if method is None:
        method = APPLY_METHOD
    if method == "exact":
        return _apply_chain(K, x, MPS, cutoff=cutoff, maxdim=maxdim)
    if method == "densitymatrix":
        return _apply_densitymatrix(K, x, cutoff=cutoff, maxdim=maxdim)
    raise ValueError("applyMPO: unknown method %r (have: exact, densitymatrix)"
                     % (method,))


def nmultMPO(A, B, cutoff=0.0, maxdim=None):
    """MPO composition. Callers are responsible for the same priming
    convention real ITensor's nmultMPO requires (see chain_session.h's
    mult_mpo(): call as nmultMPO(A, prime(B)), then mapPrime(result,2,1))
    -- this function only performs the contraction+compression, not the
    priming juggling, matching how mo_terms.h/chain_session.h keep that
    logic at the call site rather than inside the primitive."""
    return _apply_chain(A, B, MPO, cutoff=cutoff, maxdim=maxdim)


def randomMPS(sites, m, rng=None):
    """A random MPS at (at most) bond dimension m, immediately
    canonicalized (bond dimension capped by the smaller Hilbert space on
    either side of each cut, matching how an MPS's bond dimension is
    physically meaningless beyond that) and normalized. Mirrors
    mpscpp3/chain_session.h's default_mps(). `rng`, a numpy Generator,
    draws the entries from a caller's own stream; None keeps the global
    np.random stream every existing caller uses."""
    n = sites.length()
    dims = [sites.dim(i) for i in range(1, n + 1)]

    left_cum = 1
    left_bounds = []
    for d in dims[:-1]:
        left_cum = min(left_cum * d, m)
        left_bounds.append(left_cum)
    right_cum = 1
    bond_dims = list(left_bounds)
    for k in range(n - 2, -1, -1):
        right_cum = min(right_cum * dims[k + 1], m)
        bond_dims[k] = min(bond_dims[k], right_cum)

    links = [Index(d, tags="Link,l={}".format(k + 1)) for k, d in enumerate(bond_dims)]
    tensors = []
    for i in range(n):
        phys = sites.si(i + 1)
        inds = []
        if i > 0:
            inds.append(links[i - 1])
        inds.append(phys)
        if i < n - 1:
            inds.append(links[i])
        shape = tuple(ind.dim for ind in inds)
        if rng is None:
            arr = np.random.randn(*shape) + 1j * np.random.randn(*shape)
        else:
            arr = rng.standard_normal(shape) + 1j * rng.standard_normal(shape)
        tensors.append(ITensor(tuple(inds), arr))

    mps = MPS(tensors)
    mps.center = 1
    mps.position(n)  # lossless left-to-right sweep: exact canonicalization
    mps.normalize()
    return mps
