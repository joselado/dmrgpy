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
        arr = np.trace(T.array, axis1=ax_in, axis2=ax_out)
        piece = ITensor(tuple(remaining), arr)
        env = piece if env is None else env * piece
    return env.scalar()


def _bond_dims(chain):
    """{k: dim of the Link between site k and k+1} for k=1..N-1."""
    n = chain.length()
    return {k: _link_at(chain, k, k + 1).dim for k in range(1, n)}


def sum(A, B, cutoff=0.0, maxdim=None):
    """Direct sum of two MPS (or two MPO): the exact, standard MPS/MPO
    addition construction (concatenation at the boundaries, block-diagonal
    in the link spaces at interior sites), followed by one truncating
    left-to-right SVD sweep down to cutoff/maxdim."""
    n = A.length()
    if B.length() != n:
        raise ValueError("sum: mismatched chain length {} vs {}".format(n, B.length()))

    new_links = {}
    for k in range(1, n):
        la = _link_at(A, k, k + 1)
        lb = _link_at(B, k, k + 1)
        new_links[k] = Index(la.dim + lb.dim, tags="Link,l={}".format(k))

    tensors = []
    for i in range(1, n + 1):
        Ta, Tb = A.A(i), B.A(i)
        phys = tuple(ind for ind in Ta.inds if ind.hastags("Site"))
        phys_b = set(ind for ind in Tb.inds if ind.hastags("Site"))
        if set(phys) != phys_b:
            raise ValueError("sum: different index structure at site {}".format(i))

        la, ra = _link_at(A, i, i - 1), _link_at(A, i, i + 1)
        lb, rb = _link_at(B, i, i - 1), _link_at(B, i, i + 1)

        order_a = ([la] if la else []) + list(phys) + ([ra] if ra else [])
        order_b = ([lb] if lb else []) + list(phys) + ([rb] if rb else [])
        arr_a = Ta.transpose_to(order_a)
        arr_b = Tb.transpose_to(order_b)

        new_left = new_links.get(i - 1)
        new_right = new_links.get(i)
        phys_shape = tuple(p.dim for p in phys)
        shape = (((new_left.dim,) if new_left else ())
                 + phys_shape
                 + ((new_right.dim,) if new_right else ()))
        combined = np.zeros(shape, dtype=complex)

        def place(arr, left_off, left_dim, right_off, right_dim):
            idx = []
            if new_left:
                idx.append(slice(left_off, left_off + left_dim))
            idx += [slice(None)] * len(phys)
            if new_right:
                idx.append(slice(right_off, right_off + right_dim))
            combined[tuple(idx)] = arr

        place(arr_a, 0, la.dim if la else 0, 0, ra.dim if ra else 0)
        place(arr_b, la.dim if la else 0, lb.dim if lb else 0,
              ra.dim if ra else 0, rb.dim if rb else 0)

        inds = (([new_left] if new_left else []) + list(phys)
                + ([new_right] if new_right else []))
        tensors.append(ITensor(tuple(inds), combined))

    cls = type(A)
    result = cls(tensors)
    result.center = 1
    result.position(n, cutoff=cutoff, maxdim=maxdim)
    return result


def _apply_chain(K, X, out_cls, cutoff=0.0, maxdim=None):
    """Shared implementation of applyMPO (X=MPS, out_cls=MPS) and nmultMPO
    (X=MPO, out_cls=MPO): contract K against X at every site (whatever
    physical legs match, auto-contract per ITensor.__mul__), multiplying
    together each bond's K-link and X-link into one combined Link index
    (the standard, non-variational "zip up and compress" way to apply an
    MPO -- see this module's docstring), then compress with one truncating
    sweep."""
    n = K.length()
    combined_links = {}
    for k in range(1, n):
        kl = _link_at(K, k, k + 1)
        xl = _link_at(X, k, k + 1)
        combined_links[k] = Index(kl.dim * xl.dim, tags="Link,l={}".format(k))

    tensors = []
    for i in range(1, n + 1):
        Kt, Xt = K.A(i), X.A(i)
        prod = Kt * Xt

        kL, kR = _link_at(K, i, i - 1), _link_at(K, i, i + 1)
        xL, xR = _link_at(X, i, i - 1), _link_at(X, i, i + 1)
        used_links = set(l for l in (kL, kR, xL, xR) if l is not None)
        middle = [ind for ind in prod.inds if ind not in used_links]

        order = ([kL] if kL else []) + ([xL] if xL else []) + middle \
            + ([kR] if kR else []) + ([xR] if xR else [])
        arr = prod.transpose_to(order)

        left_dim = (kL.dim * xL.dim) if kL else None
        right_dim = (kR.dim * xR.dim) if kR else None
        shape = (([left_dim] if left_dim else [])
                 + [m.dim for m in middle]
                 + ([right_dim] if right_dim else []))
        arr2 = arr.reshape(tuple(shape))

        new_left = combined_links.get(i - 1) if kL else None
        new_right = combined_links.get(i) if kR else None
        inds2 = (([new_left] if new_left else []) + middle
                 + ([new_right] if new_right else []))
        tensors.append(ITensor(tuple(inds2), arr2))

    result = out_cls(tensors)
    result.center = 1
    result.position(n, cutoff=cutoff, maxdim=maxdim)
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
