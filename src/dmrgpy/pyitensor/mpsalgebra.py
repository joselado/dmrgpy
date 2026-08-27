"""Free functions over MPS/MPO chains: sum, applyMPO, nmultMPO, inner
(=innerC -- see the note below), traceC, randomMPS. Mirrors the subset of
ITensor v3's own free-function API mpscpp3/chain_session.h calls directly.

Every one of these that grows bond dimension (sum, applyMPO, nmultMPO)
finishes with a single truncating left-to-right SVD sweep (_Chain.position)
down to the requested cutoff/maxdim -- the standard, correct way to
compress a tensor-train sum/product, though not necessarily bit-for-bit
the same intermediate bond dimensions ITensor's own automaton/variational
methods would produce. That's an intentional simplification (see
autompo.py's module docstring for the same reasoning applied to MPO
construction): dmrgpy only ever observes final numerical results bounded
by Cutoff/MaxDim, never internal bond dimensions.

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
from .svd import svd
from .tensor import ITensor, commonIndex, contract_many
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


def _bond_dims(chain):
    """{k: dim of the Link between site k and k+1} for k=1..N-1."""
    n = chain.length()
    return {k: _link_at(chain, k, k + 1).dim for k in range(1, n)}


def sum_many(chains, cutoff=0.0, maxdim=None):
    """N-way direct sum of MPS (or MPO) chains: the same exact block-
    diagonal concatenation sum() does for two operands, generalized to K
    via one K-way placement per site instead of K-1 pairwise merges each
    paying their own truncating SVD sweep.

    This matters for mpobuilder.py's to_mpo(): building a T-term
    Hamiltonian by folding in one term at a time via T-1 calls to sum()
    forces T-1 full left-to-right SVD sweeps, even though (per
    mpobuilder.py's own module docstring) that per-step compression
    demonstrably doesn't reduce bond dimension at all until a final
    bidirectional pass runs -- confirmed directly, a 14-site nearest-
    neighbor Heisenberg chain (39 terms) landed at bond dimension exactly
    39, i.e. zero compression, before that final pass brought it down to
    5. Concatenating all T operands at once is exact array placement (no
    linear algebra), so only the caller's own final compression sweep(s)
    need to touch SVD at all -- this turns to_mpo()'s O(T) SVD sweeps
    into O(1)."""
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
    (X=MPO, out_cls=MPO): a single left-to-right "zip-up" sweep that
    contracts K against X one site at a time (whatever physical legs
    match, auto-contract per ITensor.__mul__) and immediately SVD-
    compresses each cut down to cutoff/maxdim before moving on, carrying
    the (already-truncated) remainder -- `leftover`, with legs (new bond,
    K's own right link, X's own right link) -- forward into the next
    site's contraction.

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
    # Right-canonicalize both inputs before the left-to-right sweep. This is
    # not a nicety, it is what makes the per-cut truncation *mean* anything:
    # at cut i the left side is already orthonormal (it is the U of the
    # previous SVD), so the singular values of `piece` are the true Schmidt
    # values of the exact product only if the not-yet-contracted right side
    # is orthonormal too. With both K and X right-orthogonal from site i+1
    # on, the fused right environment satisfies
    # sum (K†K) (x) (X†X) = I (x) I, so it is exactly orthonormal and the
    # truncation is optimal. Without it the singular values are weighted by
    # whatever norm the untouched right part happens to carry, and the cutoff
    # discards the wrong components -- confirmed directly: applying a
    # Heisenberg (H - E0 - omega) MPO to an MPS on a 10-site chain, with the
    # bond dimension (32) far *below* maxdim (60) so maxdim never even bound,
    # <Hb|Hb> and <b|H^2 b> disagreed by 0.86% (0.0018381 vs 0.0018222) where
    # the compiled ITensor backends agree to 1e-15. That ~1% error per
    # application is what left cvm.py's conjugate gradient stalling at a
    # residual of ~7e-2 on itensor_version="python" (against ~9e-6 on
    # itensor_version 2 and 3) and returning a visibly wrong spectrum.
    # This is the canonicalization step the "zip-up" algorithm assumes
    # (Stoudenmire & White, New J. Phys. 12, 055026 (2010), §3.2).
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
        U, S, V, spec = svd(piece, left_inds, cutoff=cutoff, maxdim=maxdim)
        tensors.append(U)
        leftover = S * V

    result = out_cls(tensors)
    result.center = n
    return result


def applyMPO(K, x, x0=None, cutoff=0.0, maxdim=None):
    """K (MPO) applied to x (MPS). The physical leg of the result is left
    *primed* (matching K's own "out" leg) -- callers noPrime() it back to a
    plain ket afterward, exactly like every apply_mpo() call site in
    mpscpp3/chain_session.h does.

    `x0` (an initial guess, used by real ITensor to seed an iterative
    variational "Fit" solve) is accepted for call-signature compatibility
    with chain_session.h's two applyMPO() overloads but is otherwise
    unused: this always does the same direct contract-and-compress: not
    the fastest possible method, but exact up to cutoff/maxdim regardless
    of x0's value, so ignoring it is correctness-preserving."""
    return _apply_chain(K, x, MPS, cutoff=cutoff, maxdim=maxdim)


def nmultMPO(A, B, cutoff=0.0, maxdim=None):
    """MPO composition. Callers are responsible for the same priming
    convention real ITensor's nmultMPO requires (see chain_session.h's
    mult_mpo(): call as nmultMPO(A, prime(B)), then mapPrime(result,2,1))
    -- this function only performs the contraction+compression, not the
    priming juggling, matching how mo_terms.h/chain_session.h keep that
    logic at the call site rather than inside the primitive."""
    return _apply_chain(A, B, MPO, cutoff=cutoff, maxdim=maxdim)


def randomMPS(sites, m):
    """A random MPS at (at most) bond dimension m, immediately
    canonicalized (bond dimension capped by the smaller Hilbert space on
    either side of each cut, matching how an MPS's bond dimension is
    physically meaningless beyond that) and normalized. Mirrors
    mpscpp3/chain_session.h's default_mps()."""
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
        arr = np.random.randn(*shape) + 1j * np.random.randn(*shape)
        tensors.append(ITensor(tuple(inds), arr))

    mps = MPS(tensors)
    mps.center = 1
    mps.position(n)  # lossless left-to-right sweep: exact canonicalization
    mps.normalize()
    return mps
