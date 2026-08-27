"""Two-site and one-site TDVP real-time evolution.

One time step of size dt is two half-sweeps (mirrors mpscpp3's own
tdvp_step()/TDVP/README.md convention: a left-to-right pass evolving by
dt/2, then a right-to-left pass evolving by another dt/2). At each bond in
a half-sweep (say the left-to-right one, at bond (i,i+1)):

  1. Forward-evolve the 2-site tensor theta=A(i)*A(i+1) by tau=dt/2 under
     the local 2-site effective Hamiltonian (dmrg.py's two_site_heff,
     shared with the ground-state solver -- same operator, just
     exponentiated via a short Krylov subspace instead of diagonalized).
  2. SVD-split the evolved theta into U (kept, left-orthogonal) and C=S*V.
     Set A(i) = U.
  3. If there's a next bond to process (i < N-1): *backward*-evolve the
     whole of C (bond *and* physical leg together -- a proper one-site
     tensor, dmrg.py's one_site_heff, not just the bare singular values)
     by tau, using the environment through site i (freshly rebuilt with
     the just-updated U) on one side and the *original, not-yet-touched*
     environment through site i+2 onward on the other -- crucially *not*
     an environment through site i+1, whose basis wouldn't match C's own
     bond leg (a fresh SVD-truncated one). Set A(i+1) to the result.
  4. At the sweep's last bond, skip the backward step -- C is simply
     written directly to A(i+1).

This backward correction is the whole reason two-site TDVP isn't simply
"apply the forward step at every bond and stop": it is what makes each
half-sweep equivalent (up to the truncation already being enforced by SVD)
to evolving the full N-site state by tau under all of H at once, rather
than under a sequence of independent 2-site pieces of H -- see Haegeman
et al., "Unifying time evolution and optimization with matrix product
states".

Real-time evolution needs tau purely imaginary (coeff = -i*tau forward,
+i*tau backward) -- see _lanczos_expm_multiply's docstring for the
Krylov-exponentiation method itself.

As in dmrg.py, every environment here is a <psi|...|psi> self-overlap and
built incrementally while the sweep rewrites psi's own tensors, so it
reuses dmrg.py's _extend_left/_extend_right (and their _relabel_bra_local
fix for the bra/ket link-identity collision). One extra wrinkle specific
to this module: right after SVD-splitting at bond i, site i's tensor has
a *new* link (freshly minted by svd()) that site i+1 doesn't share yet --
_extend_left/_extend_right identify a site's link by looking at its
neighbor, so each half-sweep writes a consistent (if not yet
backward-evolved) placeholder to the neighbor *before* extending the
environment past the just-updated site.

One-site TDVP (tdvp_step's num_center=1) follows the identical two-half-
sweep structure, one tensor rank lower throughout: forward-evolve psi.A(i)
alone (dmrg.py's one_site_heff, not two_site_heff -- no SVD truncation,
since one-site TDVP must conserve bond dimension exactly), split it
*losslessly* into a left/right-orthogonal site tensor and a bond tensor C,
backward-evolve C (dmrg.py's zero_site_heff -- one_site_heff's own
"backward correction" piece, one rank lower again, since a bond carries no
physical leg), and absorb it into the next site. Since one-site TDVP alone
never grows bond dimension, it needs pairing with gse.py's
global_subspace_expand() (the Yang-White Krylov global-subspace-expansion
scheme, arXiv:2005.06104/Phys. Rev. B 102, 094315) beforehand, mirroring
mpscpp3/chain_session.h's tdvp_step()/global_subspace_expand() and
TDVP/tdvp.h+TDVP/basisextension.h's own NumCenter=1 + addBasis() pairing.
"""

import numpy as np
from scipy.linalg import eigh_tridiagonal

from . import backend as bk
from .dmrg import (_all_left_environments, _all_right_environments,
                    _block_reorthogonalize, _extend_left, _extend_right,
                    _lanczos_residual, one_site_heff, two_site_heff,
                    zero_site_heff)
from .mpsalgebra import _link_at
from .svd import qr_split, svd
from .tensor import ITensor


# Convergence target for the Krylov exponentiation below, matching the
# compiled backend's own hard-wired value: ITensor's applyExp()
# (iterativesolvers.h) defaults its "ErrGoal" to 1e-10, and
# mpscpp3/chain_session.h's tdvp_step() never overrides it. Not exposed as
# a knob for the same reason `niter` isn't on that side.
_KRYLOV_ERRGOAL = 1e-10


def _lanczos_expm_multiply(matvec, v0, coeff, niter=30, tol=1e-12,
                            errgoal=_KRYLOV_ERRGOAL):
    """expm(coeff * A) @ v0 for a Hermitian linear operator A (given as a
    matvec function) via a Krylov (Lanczos) subspace of dimension up to
    `niter`: build an orthonormal basis {v0, A v0, A^2 v0, ...} via the
    Lanczos recursion (with full reorthogonalization -- the subspace is
    small enough, niter ~ tens, that this costs nothing relative to the
    matvec itself and buys real numerical stability), project A onto it as
    a small tridiagonal matrix, exponentiate it, and map back. Standard
    Krylov-propagator technique for real-time quantum dynamics; this is
    what tdvp_step()'s "niter=50 bounds the Lanczos iterations" comment in
    mpscpp3/chain_session.h refers to on the compiled-backend side.

    `niter` is an upper *bound*, not the iteration count: the recursion
    stops as soon as the Krylov approximation has converged, measured by
    Saad's a posteriori residual estimate (Y. Saad, SIAM J. Numer. Anal.
    29, 209 (1992), Sec. 5.1)

        err_k ~ ||v0|| * beta_k * |e_k^T exp(coeff T_k) e_1|,

    i.e. the size of the coupling from the last basis vector actually kept
    into the one the recursion would build next -- the same quantity the
    compiled backend's applyExp() (ITensor's iterativesolvers.h) tests via
    its Expokit-style extended-T-matrix correction. `c` (the first column of
    exp(coeff*T_k)) is needed to assemble the result anyway, so this test
    costs one extra eigh_tridiagonal of the k x k projected matrix per
    iteration -- O(k^2), against an O(D^3) MPS-level matvec -- and nothing
    else. Without it this function ran the full `niter` matvecs on *every*
    local update of *every* bond of *every* TDVP step, whether or not the
    subspace had converged 40 iterations earlier: with the callers' niter=50,
    a submode="TD" dynamical correlator on an L=30 Heisenberg chain spent
    1082 s where ~9 matvecs per call suffice.

    Reorthogonalization is done as two BLAS matrix-vector products against
    the (preallocated, incrementally filled) basis built so far, not a
    per-vector Python loop -- confirmed via cProfile that this function's
    *own* bookkeeping (independent of matvec cost) became the dominant
    remaining cost of a dynamical-METTS run after kernels.py's matvec-
    planning fix (see that module's docstring), and the O(niter) Python-
    level np.vdot calls per Lanczos step (O(niter^2) total per call here)
    were the biggest piece of it. The projected Krylov matrix is real
    tridiagonal by construction (alpha=vdot(...).real, beta=norm(...) are
    both always real), so it's exponentiated via eigh_tridiagonal's direct
    eigendecomposition rather than a generic dense expm() Pade computation
    -- cheaper, and only the first column of the exponentiated matrix is
    ever used, so the full matrix is never assembled.

    The Krylov vectors live wherever `v0` does (on the device, for a
    device backend -- see backend.py): each is chi^2 d^2 long and there are
    up to `niter` of them, so moving them is the one thing this function
    must never do. Exactly two numbers per iteration come back to the host,
    alpha and beta, because the projected tridiagonal matrix, its
    eigendecomposition and Saad's residual test are all k x k with k under
    ~10 -- host work by design, the same split dmrg.py's
    _lanczos_ground_state uses.

    Those two numbers per iteration are cheap to *move* and expensive to
    *wait for*. JAX dispatch is asynchronous: the host runs ahead
    enqueueing kernels while the device works, which is the only thing
    that hides the ~0.35 ms per-call dispatch floor this port fights
    everywhere else. `float(bk.to_host(...))` is a synchronization point,
    so reading alpha and beta every iteration drains that queue twice per
    Lanczos step -- at n=30 a single two-site TDVP step is 2 half-sweeps x
    29 bonds x ~9 iterations, i.e. ~1000 full pipeline stalls per time
    step, and a submode="TDZ" run is hundreds of time steps of exactly
    that. So on a device the recursion keeps alpha/beta as 0-d device
    arrays (`w - alpha*q` and `w/beta` need no host value) and reads a
    whole *block* of them home in one transfer -- see
    _lanczos_expm_device below, which speculates past the stopping point
    and then rolls back to the exact same k the host loop would have
    chosen, so the two paths still return the same vector."""
    _xp = bk.xp()
    beta0 = float(bk.to_host(_xp.linalg.norm(v0)))
    if beta0 == 0:
        return v0.copy()

    m = min(niter, v0.size)
    basis = _KrylovBasis(v0 / beta0, m)
    if krylov_defer_sync():
        return _lanczos_expm_device(matvec, basis, beta0, coeff, m, tol,
                                    errgoal)
    alphas = []
    betas = []

    def exp_first_column(k):
        """First column of exp(coeff * T_k), T_k the k x k projected
        (real, symmetric tridiagonal) matrix built so far. Shared with the
        device path below, which evaluates the identical quantity from a
        block of coefficients rather than one at a time."""
        return _exp_first_column(alphas, betas, coeff, k)

    beta = 0.0
    for it in range(m):
        q = basis.column(it)
        w = matvec(q)
        alpha = float(bk.to_host(_xp.vdot(q, w).real))
        alphas.append(alpha)
        if it > 0:
            w = _lanczos_residual(w, alpha, q, beta, basis.column(it - 1))
        else:
            w = w - alpha * q
        Qk = basis.matrix()  # every basis vector built so far
        w = _block_reorthogonalize(Qk, w)
        beta = float(bk.to_host(_xp.linalg.norm(w)))

        k = it + 1
        exp_col0 = exp_first_column(k)
        # beta < tol is Lanczos-sequence exhaustion (the Krylov space is
        # A-invariant and the approximation is exact); it == m-1 exhausts
        # the caller's iteration budget; otherwise stop once Saad's
        # residual estimate is below errgoal.
        if (beta < tol or it == m - 1
                or beta0 * beta * abs(exp_col0[-1]) < errgoal):
            # exp_col0 is k numbers computed on the host; this is the only
            # host->device move in the whole routine, and it is O(k).
            return beta0 * (Qk @ bk.asarray(exp_col0))

        betas.append(beta)
        basis.append(w / beta)

    raise AssertionError("unreachable: the it == m-1 branch always returns")


# How many further Lanczos iterations the device path speculates through
# before bringing a block of (alpha, beta) home. Every iteration past the
# true stopping point is a wasted matvec, and every block boundary is a
# pipeline stall, so this trades one against the other; 4 is small enough
# that the waste stays bounded even when the hint below is useless.
_KRYLOV_CHECK_BLOCK = 4

# Where to place the *first* checkpoint: the k at which the previous call
# converged. Krylov convergence here is set by the local effective
# Hamiltonian's spectral width times |coeff|, which barely moves from bond
# to bond or from step to step, so the previous k is an excellent
# predictor and the common case becomes "check once, at exactly the right
# place, waste nothing". This is a pure *scheduling* hint: it decides how
# far the recursion runs ahead, never where it stops (that is
# _first_converged_k's job, below), so a stale or wildly wrong value costs
# a few matvecs and cannot change the returned vector. Hence a plain
# process-global rather than per-chain state.
_KRYLOV_K_HINT = [1]

# Whether the Krylov exponentiator uses the deferred-synchronization path
# below. None (the default) means "on exactly when the arrays are on a
# device", which is the only case where it can help: on the host
# bk.to_host() is a no-op and speculating past the stopping point would be
# pure waste. True/False force it either way, which is what
# benchmarks/gpu/tdz_bench.py uses to attribute a device speedup to this
# change rather than to the others in the same run -- the same reason
# backend.set_jit() can be forced.
_KRYLOV_DEFER_SYNC = [None]


def set_krylov_defer_sync(mode=None):
    """Force (True) or disable (False) the deferred-synchronization Krylov
    path; None restores the automatic "on iff on a device" default. Both
    paths return the same vector -- see _lanczos_expm_device -- so this is
    a performance knob only."""
    if mode not in (True, False, None):
        raise ValueError("set_krylov_defer_sync: mode must be True, False "
                         "or None, got %r" % (mode,))
    _KRYLOV_DEFER_SYNC[0] = mode
    return mode


def krylov_defer_sync():
    """Whether the deferred-synchronization Krylov path is currently in
    use."""
    if _KRYLOV_DEFER_SYNC[0] is None:
        return bk.is_device()
    return bool(_KRYLOV_DEFER_SYNC[0])


def _exp_first_column(alphas, betas, coeff, k):
    """First column of exp(coeff * T_k), where T_k is the k x k real
    symmetric tridiagonal matrix with diagonal alphas[:k] and
    off-diagonal betas[:k-1]. Host-side by design: k is under ~10, and
    only this first column is ever used, so the full matrix exponential
    is never assembled (see _lanczos_expm_multiply's docstring)."""
    if k == 1:
        return np.array([np.exp(coeff * alphas[0])], dtype=complex)
    evals, evecs = _eigh_tridiagonal_robust(np.array(alphas[:k]),
                                             np.array(betas[:k - 1]))
    return evecs @ (np.exp(coeff * evals) * evecs[0, :])


def _first_converged_k(alphas, betas, beta0, coeff, m, tol, errgoal,
                        scan_from=1):
    """The smallest k in [scan_from, len(alphas)] at which the host loop
    in _lanczos_expm_multiply would have stopped, or None if it would
    have run past the coefficients available so far.

    This is deliberately a transcription of that loop's own three-way
    stopping test rather than an approximation of it -- beta < tol
    (Lanczos exhaustion), k == m (the caller's iteration budget), and
    Saad's a posteriori residual estimate -- because the device path's
    whole claim is that speculating ahead changes only *when* the test is
    evaluated, never *what* it decides. Returns the k, and the caller
    discards every basis vector past it."""
    for k in range(max(scan_from, 1), len(alphas) + 1):
        beta = betas[k - 1]
        if beta < tol or k == m:
            return k
        if beta0 * beta * abs(_exp_first_column(alphas, betas, coeff, k)[-1]) < errgoal:
            return k
    return None


def _lanczos_expm_device(matvec, basis, beta0, coeff, m, tol, errgoal):
    """_lanczos_expm_multiply's device path: the identical recursion with
    the per-iteration host synchronization removed.

    Two changes, and only the second one is visible from outside:

    1. alpha and beta stay 0-d device arrays. Nothing in the recursion
       needs their values on the host -- `w - alpha*q - beta*q_prev` and
       `w / beta` are array operations either way -- so keeping them where
       they are produced lets JAX's asynchronous dispatch keep running
       ahead instead of stalling twice per iteration.
    2. The stopping test therefore cannot be evaluated per iteration. So
       the recursion *speculates*: it runs to the next checkpoint, brings
       that whole block of coefficients home in one transfer, and asks
       _first_converged_k where the host loop would have stopped. If that
       point is inside the block, the extra basis vectors built past it
       are simply dropped -- they cost matvecs, they do not affect the
       answer, which is assembled from exactly the first k columns with
       exactly the same k the host path would have used.

    The one numerical wrinkle speculation introduces is that the
    recursion can step *past* a beta of ~0 (Lanczos exhaustion), where
    the host loop would already have returned. Dividing by it would put
    inf/nan into the basis and, through the reorthogonalization, into
    every later column. Clamping the divisor to 1 below instead makes
    those speculative columns exactly zero, which is harmless: they are
    orthogonal to everything, they are discarded by the rollback, and
    they cannot poison the columns that are kept."""
    _xp = bk.xp()
    alphas_d, betas_d = [], []      # 0-d device arrays, never transferred
    alphas, betas = [], []          # their host mirrors, one block at a time
    home = 0                        # coefficients already brought home
    next_check = min(max(_KRYLOV_K_HINT[0], 1), m)

    for it in range(m):
        q = basis.column(it)
        w = matvec(q)
        alpha = _xp.vdot(q, w).real
        alphas_d.append(alpha)
        if it > 0:
            w = _lanczos_residual(w, alpha, q, betas_d[it - 1],
                                  basis.column(it - 1))
        else:
            w = w - alpha * q
        Qk = basis.matrix()         # exactly k = it+1 columns
        w = _block_reorthogonalize(Qk, w)
        beta = _xp.linalg.norm(w)
        betas_d.append(beta)

        k = it + 1
        if k >= next_check or k == m:
            # One transfer for the whole block: the 2*(k-home) new
            # coefficients as a single stacked array, so a checkpoint
            # costs one synchronization rather than 2 per iteration.
            nnew = k - home
            block = bk.to_host(_xp.stack(alphas_d[home:] + betas_d[home:])).real
            alphas.extend(float(x) for x in block[:nnew])
            betas.extend(float(x) for x in block[nnew:])
            stop_k = _first_converged_k(alphas, betas, beta0, coeff, m, tol,
                                        errgoal, scan_from=home + 1)
            home = k
            if stop_k is not None:
                _KRYLOV_K_HINT[0] = stop_k
                exp_col0 = _exp_first_column(alphas, betas, coeff, stop_k)
                return beta0 * (Qk[:, :stop_k] @ bk.asarray(exp_col0))
            next_check = min(k + _KRYLOV_CHECK_BLOCK, m)

        # Clamped so a speculative step past Lanczos exhaustion produces a
        # zero column rather than nan; see this function's docstring.
        basis.append(w / _xp.where(beta > tol, beta, 1.0))

    raise AssertionError("unreachable: the k == m checkpoint always returns")


class _KrylovBasis:
    """The growing orthonormal Krylov basis of _lanczos_expm_multiply,
    stored the way the active array backend wants it stored.

    On the host it is one preallocated (n, m) buffer whose columns are
    filled in place and handed out as views -- what this function has
    always done, and what keeps the reorthogonalization two BLAS calls
    against one contiguous block.

    On a device that same buffer would be the worst available layout: JAX
    arrays are immutable, so every column assignment would copy the entire
    (n, m) block (n = chi^2 d^2, m up to 50 as the callers set it), i.e.
    O(n m) of device traffic per iteration to write O(n) of new data.
    Keeping the columns in a list and stacking the k built so far costs
    O(n k) with k typically under 10, and -- the part that actually
    matters on a device -- it is one dispatch rather than one per column.
    The stack is cached between calls to matrix() and invalidated on
    append(), so the two uses per iteration share it."""

    __slots__ = ("_on_device", "_buf", "_cols", "_k", "_stacked")

    def __init__(self, q0, m):
        self._on_device = bk.is_device()
        self._k = 1
        self._stacked = None
        if self._on_device:
            self._cols = [q0]
            self._buf = None
        else:
            self._cols = None
            self._buf = np.empty((q0.size, m), dtype=complex)
            self._buf[:, 0] = q0

    def append(self, q):
        if self._on_device:
            self._cols.append(q)
            self._stacked = None
        else:
            self._buf[:, self._k] = q
        self._k += 1

    def column(self, i):
        return self._cols[i] if self._on_device else self._buf[:, i]

    def matrix(self):
        """The (n, k) basis built so far."""
        if not self._on_device:
            return self._buf[:, :self._k]
        if self._stacked is None:
            self._stacked = bk.xp().column_stack(self._cols)
        return self._stacked


def _eigh_tridiagonal_robust(alphas, betas):
    """eigh_tridiagonal's default lapack_driver='stemr' can raise
    LinAlgError("stemr ... did not converge") on some tridiagonal
    matrices with (near-)degenerate eigenvalues -- confirmed directly on
    a real-time TDVP quench (20-orbital native-Hubbard chain, niter=50):
    reproducible at a fixed Lanczos step, not transient/load-dependent.
    'stebz' (bisection + inverse iteration) is slower but does not share
    this failure mode, so it's used as a fallback rather than letting a
    single degenerate Krylov step abort an entire quench."""
    try:
        return eigh_tridiagonal(alphas, betas)
    except np.linalg.LinAlgError:
        return eigh_tridiagonal(alphas, betas, lapack_driver="stebz")


def _evolve_two_site(L, Lbra, H, ket, i, R, Rbra, tau, niter):
    matvec, order_in, shape, x0 = two_site_heff(L, Lbra, H, ket, i, R, Rbra)
    evolved = _lanczos_expm_multiply(matvec, x0, -1j * tau, niter=niter)
    return ITensor(tuple(order_in), evolved.reshape(shape))


def _evolve_one_site(L, Lbra, H, ket, i, R, Rbra, tau, niter):
    matvec, order_in, shape, x0 = one_site_heff(L, Lbra, H, ket, i, R, Rbra)
    evolved = _lanczos_expm_multiply(matvec, x0, 1j * tau, niter=niter)
    return ITensor(tuple(order_in), evolved.reshape(shape))


def _evolve_one_site_forward(L, Lbra, H, ket, i, R, Rbra, tau, niter):
    """Forward (-i*tau) one-site evolution -- one-site TDVP's own local
    update, using the same one_site_heff() as _evolve_one_site() above
    (there used only for two-site TDVP's *backward* correction, +i*tau)."""
    matvec, order_in, shape, x0 = one_site_heff(L, Lbra, H, ket, i, R, Rbra)
    evolved = _lanczos_expm_multiply(matvec, x0, -1j * tau, niter=niter)
    return ITensor(tuple(order_in), evolved.reshape(shape))


def _evolve_zero_site(L, Lbra, R, Rbra, C, left_link, right_link, tau, niter):
    """Backward (+i*tau) evolution of a bond tensor -- one-site TDVP's
    counterpart to _evolve_one_site() above, one rank lower (see
    dmrg.py's zero_site_heff())."""
    matvec, order_in, shape, x0 = zero_site_heff(L, Lbra, R, Rbra, C, left_link, right_link)
    evolved = _lanczos_expm_multiply(matvec, x0, 1j * tau, niter=niter)
    return ITensor(tuple(order_in), evolved.reshape(shape))


def _half_sweep_lr(psi, H, tau, cutoff, maxdim, niter):
    n = psi.length()
    right_env = _all_right_environments(H, psi)  # sites i+1..N, ket = psi BEFORE this half-sweep
    left_env = {0: (None, None)}
    for i in range(1, n):
        L, Lbra = left_env[i - 1]
        R2, R2bra = right_env[i + 2]
        theta = _evolve_two_site(L, Lbra, H, psi, i, R2, R2bra, tau, niter)

        left_link = _link_at(psi, i, i - 1)
        s_i = next(ind for ind in psi.A(i).inds if ind.hastags("Site"))
        U, S, V, spec = svd(theta, ([left_link] if left_link else []) + [s_i],
                             cutoff=cutoff, maxdim=maxdim)
        psi.set_A(i, U)
        C = S * V
        # Write the (not yet backward-evolved) placeholder so site i+1
        # shares site i's freshly minted link before any neighbor-lookup
        # (_extend_left/one_site_heff's own _link_at calls) happens.
        psi.set_A(i + 1, C)
        left_env[i] = _extend_left(L, Lbra, H, psi, i)

        if i < n - 1:
            Lnew, Lnewbra = left_env[i]
            R2next, R2nextbra = right_env[i + 2]  # sites i+2..N, ORIGINAL/untouched
            C_evolved = _evolve_one_site(Lnew, Lnewbra, H, psi, i + 1, R2next, R2nextbra, tau, niter)
            psi.set_A(i + 1, C_evolved)
    psi.center = n


def _half_sweep_rl(psi, H, tau, cutoff, maxdim, niter):
    n = psi.length()
    left_env = _all_left_environments(H, psi)  # sites 1..i-1, ket = psi BEFORE this half-sweep
    right_env = {n + 1: (None, None)}
    for i in range(n - 1, 0, -1):
        L2, L2bra = left_env[i - 1]
        R, Rbra = right_env[i + 2]
        theta = _evolve_two_site(L2, L2bra, H, psi, i, R, Rbra, tau, niter)

        right_link = _link_at(psi, i + 1, i + 2)
        s_j = next(ind for ind in psi.A(i + 1).inds if ind.hastags("Site"))
        right_of_bond = [s_j] + ([right_link] if right_link else [])
        left_of_bond = [ind for ind in theta.inds if ind not in right_of_bond]
        U, S, V, spec = svd(theta, left_of_bond, cutoff=cutoff, maxdim=maxdim)
        psi.set_A(i + 1, V)
        C = U * S
        psi.set_A(i, C)
        right_env[i + 1] = _extend_right(R, Rbra, H, psi, i + 1)

        if i > 1:
            L2prev, L2prevbra = left_env[i - 1]  # sites 1..i-1, ORIGINAL/untouched
            Rnew, Rnewbra = right_env[i + 1]
            C_evolved = _evolve_one_site(L2prev, L2prevbra, H, psi, i, Rnew, Rnewbra, tau, niter)
            psi.set_A(i, C_evolved)
    psi.center = 1


def _half_sweep_lr_onesite(psi, H, tau, niter):
    """One-site analogue of _half_sweep_lr() above: at each site i,
    forward-evolve psi.A(i) alone (one_site_heff, not two_site_heff -- no
    truncation, since one-site TDVP must conserve bond dimension exactly),
    split it *losslessly* via QR (svd.py's qr_split -- cheaper than a full
    SVD since no truncation/singular values are needed here) into a
    left-orthogonal Q kept at site i and a bond tensor C, then
    backward-evolve C (zero_site_heff, dmrg.py) and absorb
    it into site i+1 before that site's own forward step. Mirrors
    TDVP/tdvp.h's NumCenter=1 sweep. Pair with global_subspace_expand()
    (gse.py) beforehand -- this alone never grows bond dimension."""
    n = psi.length()
    right_env = _all_right_environments(H, psi)  # sites i+1..N, BEFORE this half-sweep
    left_env = {0: (None, None)}
    for i in range(1, n + 1):
        L, Lbra = left_env[i - 1]
        Rn, Rnbra = right_env[i + 1]
        left_link = _link_at(psi, i, i - 1)
        right_link = _link_at(psi, i, i + 1)
        A_new = _evolve_one_site_forward(L, Lbra, H, psi, i, Rn, Rnbra, tau, niter)
        if i == n:
            psi.set_A(i, A_new)
            continue
        s_i = next(ind for ind in A_new.inds if ind.hastags("Site"))
        Q, C, new_link = qr_split(A_new, ([left_link] if left_link else []) + [s_i],
                                   orthonormal="left")
        psi.set_A(i, Q)
        orig_next = psi.A(i + 1)
        # Write a consistent (if not yet backward-evolved) placeholder to
        # site i+1 *before* extending the environment past site i: Q's own
        # new link isn't shared with orig_next yet, and _extend_left looks
        # up a site's link by its neighbor (mirrors _half_sweep_lr's own
        # same wrinkle, see this module's docstring -- except there C is
        # already a full site-shaped tensor from splitting a 2-site blob,
        # so it doubles directly as the placeholder; here C is bond-only,
        # so C*orig_next is used instead).
        psi.set_A(i + 1, C * orig_next)
        left_env[i] = _extend_left(L, Lbra, H, psi, i)
        Lnew, Lnewbra = left_env[i]
        C_evolved = _evolve_zero_site(Lnew, Lnewbra, Rn, Rnbra, C, new_link, right_link, tau, niter)
        psi.set_A(i + 1, C_evolved * orig_next)
    psi.center = n


def _half_sweep_rl_onesite(psi, H, tau, niter):
    """Mirror of _half_sweep_lr_onesite() above, sweeping right to left."""
    n = psi.length()
    left_env = _all_left_environments(H, psi)  # sites 1..i-1, BEFORE this half-sweep
    right_env = {n + 1: (None, None)}
    for i in range(n, 0, -1):
        L2, L2bra = left_env[i - 1]
        R, Rbra = right_env[i + 1]
        left_link = _link_at(psi, i, i - 1)
        right_link = _link_at(psi, i, i + 1)
        A_new = _evolve_one_site_forward(L2, L2bra, H, psi, i, R, Rbra, tau, niter)
        if i == 1:
            psi.set_A(i, A_new)
            continue
        s_i = next(ind for ind in A_new.inds if ind.hastags("Site"))
        right_of_bond = [s_i] + ([right_link] if right_link else [])
        left_of_bond = [ind for ind in A_new.inds if ind not in right_of_bond]
        C, V, new_link = qr_split(A_new, left_of_bond, orthonormal="right")
        psi.set_A(i, V)
        orig_prev = psi.A(i - 1)
        # Placeholder write before extending past site i -- see
        # _half_sweep_lr_onesite's matching comment.
        psi.set_A(i - 1, orig_prev * C)
        right_env[i] = _extend_right(R, Rbra, H, psi, i)
        Rnew, Rnewbra = right_env[i]
        C_evolved = _evolve_zero_site(L2, L2bra, Rnew, Rnewbra, C, left_link, new_link, tau, niter)
        psi.set_A(i - 1, orig_prev * C_evolved)
    psi.center = 1


def tdvp_step(psi, H, dt, cutoff, maxdim, niter=50, num_center=2):
    """One real-time step exp(-i*dt*H) via TDVP: a left-to-right half-sweep
    evolving by dt/2, then a right-to-left half-sweep evolving by another
    dt/2 -- mirrors mpscpp3/chain_session.h's tdvp_step(). num_center=2
    (default, matches every pre-existing caller) runs two-site TDVP,
    which grows bond dimension via each bond's SVD truncation (cutoff,
    maxdim used); num_center=1 runs one-site TDVP, which conserves bond
    dimension exactly (cutoff/maxdim unused -- pair with
    gse.global_subspace_expand() to grow it beforehand, the Yang-White
    scheme this module's own one-site path is meant to be paired with).
    Mutates psi in place."""
    tau = dt / 2.0
    if num_center == 1:
        _half_sweep_lr_onesite(psi, H, tau, niter)
        _half_sweep_rl_onesite(psi, H, tau, niter)
    else:
        _half_sweep_lr(psi, H, tau, cutoff, maxdim, niter)
        _half_sweep_rl(psi, H, tau, cutoff, maxdim, niter)
    return psi
