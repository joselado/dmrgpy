"""Two-site DMRG ground-state solver.

Standard two-site sweep: build left/right environments of <psi|H|psi>,
diagonalize the local two-site effective Hamiltonian at each bond via a
hand-rolled Lanczos ground-state solver (_lanczos_ground_state below --
the 2-site block is never densified beyond what its own matvec needs),
SVD-split the result back into two site tensors with truncation, move the
orthogonality center, and repeat sweeps per the given Sweeps schedule.

_lanczos_ground_state() replaces an earlier scipy.sparse.linalg.eigsh
(ARPACK) call. Measured directly via cProfile on a representative sweep,
ARPACK's own per-call iteration/dispatch bookkeeping
(scipy/sparse/linalg/eigen/arpack/arpack.py's iterate loop) accounted for
1.45s of tottime by itself, separate from and on top of the matvec calls
it drives -- overhead intrinsic to going through ARPACK's reverse-
communication Fortran interface for a problem this small (k=1, small
Krylov dimension), not to the eigenproblem itself. The replacement
mirrors tdvp.py's _lanczos_expm_multiply's Krylov construction (Lanczos
recursion with full reorthogonalization -- the subspace is small enough,
niter ~ tens, for this to cost nothing and buy real numerical stability),
except it tracks the lowest Ritz value's own convergence (stopping once
it stabilizes) rather than running a fixed number of steps, since unlike
expm_multiply's smooth exponential map, a ground-state eigenvalue can
stabilize well before niter steps -- exiting early both matches ARPACK's
own tol-based early exit and skips wasted matvecs.

The one subtlety worth a comment (see _relabel_bra_local): every
environment here is a <psi|...|psi>-style self-overlap, so it hits the
exact same bra/ket link-identity collision mpsalgebra.py's inner() had to
work around (_fresh_link_copy) -- except here the environment is built
*incrementally*, one site at a time, interleaved with a sweep that
actively rewrites psi's own tensors as it goes, so the fix has to be
applied per-site-as-we-go rather than once upfront.

Noise-term support (density-matrix perturbation, meant to help DMRG escape
symmetric local minima) is not implemented -- randomMPS()'s starting point
already avoids the specific trap that exists to fix (see get_sites.h's
comment on why mpscpp3 seeds from an actual random MPS rather than a
product state), and dropping it is exactly the kind of optionality
CLAUDE.md's plan allows as long as dmrgpy's own results aren't affected.
Sweeps.noise is accepted (for call-signature compatibility with
mpscpp3/chain_session.h's make_sweeps()) but otherwise unused.

dmrg_generalized() below solves the *generalized* eigenproblem
H|psi>=lambda*A|psi> (A Hermitian positive definite) by reusing the same
per-sweep machinery (_dmrg_one_sweep, factored out of dmrg()'s own loop
body) against a shifted operator H-lambda*A, self-consistently updating
lambda between sweeps -- see its own docstring for the derivation.
"""

import numpy as np

from . import backend as bk

from . import kernels
from .mpsalgebra import _link_at, inner
from .mpsalgebra import sum as mps_sum
from .svd import svd
from .tensor import ITensor, contract_many, dag
from .tensor import prime as _t_prime


def _tridiag_matrix(alphas, betas):
    """The dense symmetric tridiagonal Lanczos matrix, built with np.diag
    rather than an element-by-element Python loop."""
    a = np.asarray(alphas, dtype=float)
    if len(a) == 1:
        return a.reshape(1, 1)
    b = np.asarray(betas, dtype=float)
    return np.diag(a) + np.diag(b, 1) + np.diag(b, -1)


def _tridiag_ground_value(alphas, betas):
    """Lowest eigenVALUE only of the Lanczos tridiagonal matrix -- all the
    per-iteration convergence test needs. The eigenvector is built once, on
    the iteration that actually returns.

    `eigvalsh` rather than `eigh` because the eigenvectors were being
    computed and discarded on every iteration (5799 of them in a 30-site
    maxm=30 solve). Benchmarked at the k that actually occur here (k=3..20,
    the Lanczos subspace rarely grows past ~10): this form costs 17-39us
    against 19-90us for the original loop-built `eigh`. scipy's dedicated
    `eigvalsh_tridiagonal` was tried and is *worse* at these sizes (39-76us)
    -- its Python-level argument-validation wrapper costs more than the O(k^3)
    it saves, which only pays off for k well beyond anything Lanczos reaches
    here."""
    return float(np.linalg.eigvalsh(_tridiag_matrix(alphas, betas))[0])


def _tridiag_ground_ritz(alphas, betas):
    """Lowest eigenpair of the real symmetric tridiagonal Lanczos matrix
    built from `alphas` (diagonal) and `betas` (off-diagonal) so far."""
    w, v = np.linalg.eigh(_tridiag_matrix(alphas, betas))
    return w[0], v[:, 0]


def _lanczos_residual_impl(w, alpha, q, beta, q_prev):
    """The Lanczos three-term recurrence, w <- w - alpha q_k - beta q_{k-1},
    as one function so it fuses into a single kernel under `backend.jit`
    (three eager dispatches otherwise). alpha/beta stay traced values, not
    static ones -- specializing on their *values* would retrace every
    iteration."""
    return w - alpha * q - beta * q_prev


_lanczos_residual = bk.jit(_lanczos_residual_impl)


def _block_reorthogonalize_impl(Q, w):
    """w minus its projection onto the span of Q's columns, as two matmuls
    (one kernel under jit) instead of one dispatch per basis vector."""
    return w - Q @ (Q.conj().T @ w)


_block_reorthogonalize = bk.jit(_block_reorthogonalize_impl)


def _lanczos_ground_state(matvec, v0, niter=30, tol=1e-12):
    """Lowest eigenpair of a Hermitian linear operator (given as a matvec
    function) via Lanczos with full reorthogonalization, stopping early
    once the lowest Ritz value stabilizes to `tol` (relative). See this
    module's docstring for why this replaces scipy.sparse.linalg.eigsh
    here. Returns (eigenvalue, eigenvector) with eigenvector normalized
    (v0 need not be)."""
    # The Krylov vectors stay wherever v0 lives (on the device, for a
    # device backend). Only the tridiagonal coefficients alpha/beta come
    # back to the host -- two numbers per iteration, against O(chi^2 d^2)
    # of vector that never moves. The tridiagonal eigenproblem itself is
    # at most 30x30 and runs on the host by design.
    _xp = bk.xp()
    on_device = bk.is_device()
    beta0 = float(bk.to_host(_xp.linalg.norm(v0)))
    if beta0 == 0:
        v0 = bk.asarray(np.random.default_rng(0).standard_normal(v0.shape) + 0j)
        beta0 = float(bk.to_host(_xp.linalg.norm(v0)))
    q = v0 / beta0
    qs = [q]
    alphas = []
    betas = []
    w = matvec(q)
    alpha = float(bk.to_host(_xp.vdot(q, w).real))
    alphas.append(alpha)
    w = w - alpha * q

    prev_eval = _tridiag_ground_value(alphas, betas)
    m = min(niter, v0.size)
    for _ in range(1, m):
        beta = float(bk.to_host(_xp.linalg.norm(w)))
        if beta < tol:
            break
        betas.append(beta)
        q_new = w / beta
        qs.append(q_new)
        w = matvec(q_new)
        alpha = float(bk.to_host(_xp.vdot(q_new, w).real))
        alphas.append(alpha)
        w = _lanczos_residual(w, alpha, q_new, beta, qs[-2])
        if on_device:
            # One block projection instead of one dispatch per basis
            # vector: at the subspace sizes Lanczos reaches here (k up to
            # 30) the loop below issues 2k tiny kernels, each paying the
            # ~0.35 ms device dispatch floor -- more, at small chi, than
            # the matvec they orthogonalize against. The host keeps the
            # loop: it is modified Gram-Schmidt, marginally more stable
            # than this classical/block form, it costs nothing there, and
            # the NumPy path of this port is meant to stay bit-identical
            # to what it was before the port existed.
            w = _block_reorthogonalize(_xp.column_stack(qs[:-1]), w)
        else:
            for qk in qs[:-1]:
                # the overlap stays on the device: it is immediately
                # multiplied back into a device vector, so pulling it home
                # would be a round trip for nothing
                w = w - _xp.vdot(qk, w) * qk

        # eigenvalue only for the convergence test; the eigenvector is
        # built once, on the iteration that actually returns
        cur_eval = _tridiag_ground_value(alphas, betas)
        if abs(cur_eval - prev_eval) < tol * max(1.0, abs(cur_eval)):
            cur_eval, cur_vec = _tridiag_ground_ritz(alphas, betas)
            Q = _xp.column_stack(qs)
            return cur_eval, Q @ bk.asarray(cur_vec)
        prev_eval = cur_eval

    cur_eval, cur_vec = _tridiag_ground_ritz(alphas, betas)
    Q = _xp.column_stack(qs)
    return cur_eval, Q @ bk.asarray(cur_vec)


def _relabel_bra_local(T, chain, i, left_bra, right_bra):
    """chain.A(i), turned into one more bra-local piece of a left/right
    environment: its own Link legs are relabeled (left -> left_bra if
    given, else a fresh Index; same for right), then its physical leg is
    primed (to match H's output leg) and the tensor conjugated. Mirrors
    mpsalgebra.py's _fresh_link_copy, done one site at a time so it can be
    interleaved with a sweep that rewrites chain's own tensors as it goes."""
    old_left = _link_at(chain, i, i - 1)
    old_right = _link_at(chain, i, i + 1)
    if old_left is not None and left_bra is None:
        left_bra = old_left.sim()
    if old_right is not None and right_bra is None:
        right_bra = old_right.sim()
    new_inds = []
    for ind in T.inds:
        if old_left is not None and ind == old_left:
            new_inds.append(left_bra)
        elif old_right is not None and ind == old_right:
            new_inds.append(right_bra)
        else:
            new_inds.append(ind)
    bra_piece = dag(_t_prime(ITensor(tuple(new_inds), T.array), "Site"))
    return bra_piece, left_bra, right_bra


def _extend_left(L, left_bra, H, ket, i, bra=None):
    """One more site of the left environment <bra|H|ket>. `bra` defaults
    to `ket` (the ordinary self-overlap environments of ground-state
    DMRG); nhdmrg.py passes a *different* MPS as bra (its two-sided
    environments), which works unchanged because the fidelity-truncated
    left/right MPS share link Index objects, hitting the exact same
    bra/ket link-identity collision _relabel_bra_local exists to solve."""
    if bra is None:
        bra = ket
    T = ket.A(i)
    bra_piece, _, right_bra = _relabel_bra_local(bra.A(i), bra, i, left_bra, None)
    # contract_many(), not a left-to-right `piece = bra_piece*H.A(i)*T;
    # new_L = piece if L is None else L*piece` chain: that ordering builds
    # bra_piece*H.A(i)*T *before* L ever gets a chance to cancel away one
    # side of bra_piece's/T's own link legs, so the intermediate carries
    # both of them simultaneously -- measured directly at maxdim=60 on a
    # 14-site chain, a 92-million-element intermediate (1.5s) instead of
    # a few thousand elements (0.0006s) for the mathematically identical
    # contract_many() result. See tensor.py's contract_many() docstring.
    pieces = [p for p in (L, bra_piece, H.A(i), T) if p is not None]
    new_L = contract_many(pieces)
    return new_L, right_bra


def _extend_right(R, right_bra, H, ket, i, bra=None):
    if bra is None:
        bra = ket
    T = ket.A(i)
    bra_piece, left_bra, _ = _relabel_bra_local(bra.A(i), bra, i, None, right_bra)
    pieces = [p for p in (R, bra_piece, H.A(i), T) if p is not None]
    new_R = contract_many(pieces)
    return new_R, left_bra


def _all_right_environments(H, ket, bra=None):
    """{i: (R_tensor_or_None, dangling_bra_link_or_None)} for i = N+1..2."""
    n = ket.length()
    env = {n + 1: (None, None)}
    for i in range(n, 1, -1):
        R_next, bra_next = env[i + 1]
        env[i] = _extend_right(R_next, bra_next, H, ket, i, bra=bra)
    return env


def _all_left_environments(H, ket):
    """{i: (L_tensor_or_None, dangling_bra_link_or_None)} for i = 0..N-1."""
    n = ket.length()
    env = {0: (None, None)}
    for i in range(1, n):
        L_prev, bra_prev = env[i - 1]
        env[i] = _extend_left(L_prev, bra_prev, H, ket, i)
    return env


def two_site_heff(L, Lbra, H, ket, i, R, Rbra):
    """The 2-site effective Hamiltonian at bond (i,i+1), as a matvec over
    flat arrays shaped like ket.A(i)*ket.A(i+1) -- shared by dmrg.py's own
    Lanczos ground-state search and tdvp.py's Krylov real-time propagation,
    so both build the identical operator. Returns (matvec, order_in,
    shape, x0) where x0 is ket's own current 2-site tensor, flattened in
    `order_in`'s axis order (which the caller reshapes results back into
    to build an ITensor)."""
    Ti, Tj = ket.A(i), ket.A(i + 1)
    left_link = _link_at(ket, i, i - 1)
    right_link = _link_at(ket, i + 1, i + 2)
    s_i = next(ind for ind in Ti.inds if ind.hastags("Site"))
    s_j = next(ind for ind in Tj.inds if ind.hastags("Site"))

    order_in = ([left_link] if left_link else []) + [s_i, s_j] + ([right_link] if right_link else [])
    shape = tuple(ind.dim for ind in order_in)
    x0 = (Ti * Tj).transpose_to(order_in).reshape(-1)

    s_i_out, s_j_out = s_i.prime(1), s_j.prime(1)
    order_out = ([Lbra] if Lbra else []) + [s_i_out, s_j_out] + ([Rbra] if Rbra else [])
    H_i, H_j = H.A(i), H.A(i + 1)
    pieces = [p for p in (L, H_i, H_j, R) if p is not None]
    matvec = kernels.make_matvec(pieces, order_in, shape, order_out)

    return matvec, order_in, shape, x0


def one_site_heff(L, Lbra, H, ket, i, R, Rbra):
    """The 1-site effective Hamiltonian at site i alone -- the "backward
    evolution" piece that is specifically what makes two-site TDVP
    different from two independent single-bond time-evolution steps (see
    tdvp.py's module docstring). Acts on ket.A(i) itself (bond * physical
    leg together, e.g. S*V from a neighboring bond's SVD split, temporarily
    written into ket.A(i) by the caller so link identities resolve
    consistently -- see tdvp.py's half-sweep functions). Same shape of
    return value as two_site_heff()."""
    T = ket.A(i)
    left_link = _link_at(ket, i, i - 1)
    right_link = _link_at(ket, i, i + 1)
    s_i = next(ind for ind in T.inds if ind.hastags("Site"))

    order_in = ([left_link] if left_link else []) + [s_i] + ([right_link] if right_link else [])
    shape = tuple(ind.dim for ind in order_in)
    x0 = T.transpose_to(order_in).reshape(-1)

    s_i_out = s_i.prime(1)
    order_out = ([Lbra] if Lbra else []) + [s_i_out] + ([Rbra] if Rbra else [])
    H_i = H.A(i)
    pieces = [p for p in (L, H_i, R) if p is not None]
    matvec = kernels.make_matvec(pieces, order_in, shape, order_out)

    return matvec, order_in, shape, x0


def zero_site_heff(L, Lbra, R, Rbra, C, left_link, right_link):
    """The 0-site ("bond") effective Hamiltonian sitting between two
    already-extended environments L (through some site i) and R (through
    site i+1 onward), with no local MPO tensor of its own in between since
    a bond carries no physical index -- the backward-evolution piece that
    makes one-site TDVP's forward site-update sweep equivalent to
    evolving the whole state by tau under all of H at once (Haegeman et
    al., "Unifying time evolution and optimization with matrix product
    states"), exactly one_site_heff's role in two-site TDVP but one rank
    lower. C is the current bond tensor (indices left_link, right_link
    only); left_link must be the index shared with L's own dangling ket
    leg, right_link the one shared with R's. Same shape of return value
    as one_site_heff/two_site_heff."""
    order_in = ([left_link] if left_link else []) + ([right_link] if right_link else [])
    shape = tuple(ind.dim for ind in order_in)
    x0 = C.transpose_to(order_in).reshape(-1)

    order_out = ([Lbra] if Lbra else []) + ([Rbra] if Rbra else [])
    pieces = [p for p in (L, R) if p is not None]
    matvec = kernels.make_matvec(pieces, order_in, shape, order_out)

    return matvec, order_in, shape, x0


def _local_ground_state(L, Lbra, R, Rbra, H, ket, i, niter):
    """Diagonalize the 2-site effective Hamiltonian at bond (i,i+1) for its
    lowest eigenpair. Returns (energy, theta_ITensor)."""
    matvec, order_in, shape, x0 = two_site_heff(L, Lbra, H, ket, i, R, Rbra)

    dim = x0.size
    if dim <= 3:
        # too small a space for Lanczos to be meaningful; diagonalize directly.
        # dim <= 3, so this whole matrix is at most 3x3: build it with the
        # active namespace (matvec returns device arrays for a device
        # backend) and diagonalize it there rather than transferring.
        _xp = bk.xp()
        basis = bk.eye(dim)
        Hmat = _xp.column_stack([matvec(basis[:, k]) for k in range(dim)])
        w, v = _xp.linalg.eigh((Hmat + Hmat.conj().T) / 2)
        eval0, evec0 = w[0], v[:, 0]
    else:
        # The Sweeps schedule's own niter (e.g. ITensor's usual default of
        # 2) is deliberately small -- it relies on DMRG's many sweeps to
        # refine the state rather than fully converging each local
        # diagonalization -- but that only works once the sweep is already
        # close to converged; the OLD eigsh call floored this at
        # max(niter,200) (via maxiter=) for exactly this reason, and
        # dropping that floor when eigsh was replaced was a real
        # regression, confirmed directly: without it, larger local Hilbert
        # spaces (e.g. spin-3/2) and the excited-states penalty method
        # stopped converging to the true ground/target state (up to 4.7e-4
        # energy error on a case that should match ED to ~1e-14, and a
        # first-excited-state gap off by 2.7x). _lanczos_ground_state's own
        # early-stop-on-stable-Ritz-value logic means this floor costs
        # little in practice -- most bonds still exit in far fewer than
        # 200 iterations once bond dimension saturates.
        eval0, evec0 = _lanczos_ground_state(matvec, x0, niter=max(niter, 200))

    theta = ITensor(tuple(order_in), evec0.reshape(shape))
    return float(eval0.real), theta


def _apply_local_update(ket, i, theta, cutoff, maxdim, direction):
    """SVD-split theta (spanning sites i,i+1) back into ket.A(i)/ket.A(i+1)
    with truncation, moving the orthogonality center in `direction`
    ('right': keep U left-orthogonal at i, absorb S*V into i+1; 'left':
    mirror)."""
    left_link = _link_at(ket, i, i - 1)
    s_i = next(ind for ind in ket.A(i).inds if ind.hastags("Site"))
    left_inds = ([left_link] if left_link else []) + [s_i]
    U, S, V, spec = svd(theta, left_inds, cutoff=cutoff, maxdim=maxdim)
    if direction == "right":
        ket.set_A(i, U)
        ket.set_A(i + 1, S * V)
        ket.center = i + 1
    else:
        ket.set_A(i, U * S)
        ket.set_A(i + 1, V)
        ket.center = i


def _dmrg_one_sweep(psi, H, maxdim, cutoff, niter):
    """One full right-then-left two-site ground-state sweep against a
    fixed H, mutating psi in place. Returns the last local bond's
    eigenvalue (dmrg()'s own returned-`energy` convention -- see its
    docstring). Factored out of dmrg()'s own per-sweep loop body so
    dmrg_generalized() below can reuse the identical sweep machinery
    while rebuilding H itself between sweeps."""
    n = psi.length()
    energy = None

    right_env = _all_right_environments(H, psi)
    left_env = {0: (None, None)}
    for i in range(1, n):
        L, Lbra = left_env[i - 1]
        R, Rbra = right_env[i + 2]
        energy, theta = _local_ground_state(L, Lbra, R, Rbra, H, psi, i, niter)
        _apply_local_update(psi, i, theta, cutoff, maxdim, "right")
        left_env[i] = _extend_left(L, Lbra, H, psi, i)

    right_env = {n + 1: (None, None)}
    for i in range(n - 1, 0, -1):
        L, Lbra = left_env[i - 1]
        R, Rbra = right_env[i + 2]
        energy, theta = _local_ground_state(L, Lbra, R, Rbra, H, psi, i, niter)
        _apply_local_update(psi, i, theta, cutoff, maxdim, "left")
        right_env[i + 1] = _extend_right(R, Rbra, H, psi, i + 1)

    return energy


def dmrg(psi, H, sweeps, quiet=True):
    """Two-site ground-state DMRG. Mutates psi in place and returns the
    final energy -- mirrors ITensor's own dmrg(psi,H,sweeps,args) signature
    (minus the Args bag, which this package uses explicit kwargs for
    throughout, and minus the excited-state penalty overload -- see
    dmrg_excited() for that)."""
    energy = None
    for sweep_i in range(sweeps.nsweep):
        maxdim, cutoff, noise, niter = sweeps.at(sweep_i)
        energy = _dmrg_one_sweep(psi, H, maxdim, cutoff, niter)
        if not quiet:
            print("sweep {}: energy = {}".format(sweep_i, energy))

    return energy


def lam0_is_unset(lam0):
    """True if lam0 (a dmrg_generalized()/nhdmrg_generalized() starting
    lambda estimate) means "unset": None, or a float/complex/numpy scalar
    whose real part is NaN. Both a bare None and a NaN sentinel are
    accepted (see dmrg_generalized()'s own docstring for why NaN needs to
    be accepted at all -- the C++/pybind11 session methods have no None
    equivalent). Coerces via complex(lam0) rather than gating on a
    specific type (isinstance(..., float) or isinstance(..., complex))
    precisely because either gate misses the *other* numeric type: a
    caller passing a bare float NaN to nhdmrg_generalized() (whose own
    check used to require isinstance(..., complex)) or a complex lam0 to
    dmrg_generalized() (whose own check used to require
    isinstance(..., float)) would otherwise silently skip the auto-seed
    fallback and feed a live NaN straight into the first outer
    iteration's Heff construction -- confirmed directly, found via code
    review."""
    if lam0 is None:
        return True
    try:
        return np.isnan(complex(lam0).real)
    except (TypeError, ValueError):
        return False


def dmrg_generalized(psi, H, A, sweeps, quiet=True, lam0=None):
    """Generalized-eigenvalue ground-state DMRG: finds the smallest
    lambda solving $H|\\psi\\rangle=\\lambda A|\\psi\\rangle$ for a Hermitian
    positive-definite metric MPO A (A = the identity MPO reduces this
    exactly to plain dmrg()). Mutates psi in place and returns the final
    lambda. H itself is assumed Hermitian too (same precondition plain
    dmrg() already has -- unlike gs_energy(), which dispatches non-
    Hermitian H to NH-DMRG, nothing at this level checks this, and the
    local two-site solver silently Hermitizes its effective-Hamiltonian
    matrix before diagonalizing regardless; groundstate.py's Python-level
    wrapper is where this actually gets validated for real callers).

    Unlike Chain::gs_energy_generalized's ITensor v3 C++ port
    (mpscpp3/chain_session.h), this pure-Python engine never implements
    DMRG's noise term at all (see this module's own docstring) -- so for
    a Hamiltonian/metric pair whose self-consistent iteration has a
    symmetric local minimum, the two backends are not guaranteed to
    converge to the same lambda beyond ordinary numerical tolerance
    (the v3 port can use noise to escape such a minimum; this one
    cannot).

    The trick: minimizing $\\langle\\psi|H|\\psi\\rangle$ subject to the
    *metric* normalization $\\langle\\psi|A|\\psi\\rangle=1$ has stationarity
    condition $(H-\\lambda A)|\\psi\\rangle=\\mu|\\psi\\rangle$ for Lagrange
    multiplier lambda -- the ordinary (mu,psi) eigenproblem of the
    *plain*-normalized ($\\langle\\psi|\\psi\\rangle=1$) shifted operator
    $H-\\lambda A$, exactly what a standard two-site DMRG sweep already
    finds. At $\\mu=0$ this is precisely $H|\\psi\\rangle=\\lambda
    A|\\psi\\rangle$, so each outer iteration here (i) builds the MPO
    $H_\\mathrm{eff}=H-\\lambda A$ from the current lambda estimate, (ii)
    runs one ordinary DMRG sweep (_dmrg_one_sweep, unmodified) against
    $H_\\mathrm{eff}$, then (iii) updates lambda to the freshly-swept
    state's generalized Rayleigh quotient
    $\\langle\\psi|H|\\psi\\rangle/\\langle\\psi|A|\\psi\\rangle$ -- a
    self-consistent-field-style fixed-point iteration on lambda that also
    drives mu to 0 at convergence. One outer iteration per Sweeps
    schedule entry, so bond dimension/cutoff/niter ramp exactly as an
    ordinary dmrg() run's own schedule specifies. $H_\\mathrm{eff}$ is
    rebuilt from the *original* H, A each outer iteration (never
    accumulated), so its own bond dimension never grows across outer
    iterations -- always at most bondDim(H)+bondDim(A), summed exactly
    (cutoff=0, maxdim=None: H and A are typically already-compressed,
    modest-bond-dimension Hamiltonian-like MPOs, so an exact sum costs
    little and avoids silently truncating the shift term).

    lam0 (optional): starting lambda estimate; None or NaN both mean
    "unset" (NaN accepted too since Chain::gs_energy_generalized's own
    C++/pybind11 binding has no None equivalent and uses NaN as its own
    "unset" sentinel -- accepting both here keeps the two session-level
    APIs interchangeable). Unset defaults to psi's own current
    generalized Rayleigh quotient $\\langle\\psi|H|\\psi\\rangle/
    \\langle\\psi|A|\\psi\\rangle$, which costs two cheap inner()s and gives
    the first outer sweep a head start over starting from lambda=0.

    Raises RuntimeError if <psi|A|psi> ever collapses to ~0 mid-iteration
    (only expected if A isn't actually positive definite, or a sweep
    drives psi toward A's near-null-space) -- dividing through such a
    value would otherwise silently poison every later outer iteration
    with inf/NaN instead of failing loudly.
    """
    lam = lam0
    if lam0_is_unset(lam):
        a0 = inner(psi, A, psi).real
        h0 = inner(psi, H, psi).real
        lam = h0 / a0 if abs(a0) > 1e-14 else 0.0

    for sweep_i in range(sweeps.nsweep):
        maxdim, cutoff, noise, niter = sweeps.at(sweep_i)
        Heff = mps_sum(H, A * (-lam))
        _dmrg_one_sweep(psi, Heff, maxdim, cutoff, niter)
        a_psi = inner(psi, A, psi).real
        if not (abs(a_psi) > 1e-14):
            # written as "not (> tol)" rather than "< tol": a NaN a_psi
            # (e.g. from numerical blowup elsewhere) fails *both*
            # comparisons, so the more obvious "< tol" form would let a
            # NaN silently slip past this guard and propagate as a
            # seemingly-converged NaN lambda -- confirmed by code review.
            raise RuntimeError(
                "dmrg_generalized: <psi|A|psi> collapsed to ~0 (or NaN) "
                "mid-iteration (A may not be positive definite, or the "
                "sweep drove psi toward A's near-null-space) -- cannot "
                "form a meaningful generalized Rayleigh quotient")
        h_psi = inner(psi, H, psi).real
        lam = h_psi / a_psi
        if not quiet:
            print("sweep {}: lambda = {}".format(sweep_i, lam))

    return lam


# -- excited states: overlap-penalty method --------------------------------
#
# mirrors chain_session.h's excited_states(): dmrg(psi1,H_,wfs,sweeps,args)
# with a "Weight" arg augments the objective with weight*sum_k|wfs_k><wfs_k|,
# penalizing overlap with already-found lower states. Since there's no H
# operator mediating a projector term, its 2-site-local contribution is
# just a plain overlap: for a fixed already-found state wfs_k, precompute
# proj_k = <wfs_k(left)|psi(left)> * wfs_k's own 2-site tensor *
# <wfs_k(right)|psi(right)> once per bond (independent of the Krylov
# vector) -- since proj_k ends up expressed in exactly psi's own current
# link+physical indices (see the derivation this was checked against in
# selftest_dmrg.py), its contribution to the effective-Hamiltonian matvec
# is just weight * <proj_k|v> * proj_k, a rank-1 update computed directly
# on flat numpy arrays alongside the ordinary two_site_heff matvec.
#
# No bra/ket link-identity collision risk here (unlike the H-environments'
# _relabel_bra_local dance): wfs_k and psi are always genuinely different
# MPS objects with independently-minted link ids, so a plain dag() on
# wfs_k's tensors is safe.

def _extend_overlap_left(L, wfs_k, psi, i):
    piece = dag(wfs_k.A(i)) * psi.A(i)
    return piece if L is None else L * piece


def _extend_overlap_right(R, wfs_k, psi, i):
    piece = dag(wfs_k.A(i)) * psi.A(i)
    return piece if R is None else piece * R


def _all_overlap_right(wfs_k, psi):
    n = psi.length()
    env = {n + 1: None}
    for i in range(n, 1, -1):
        env[i] = _extend_overlap_right(env[i + 1], wfs_k, psi, i)
    return env


def _bond_projections(penalty_states, left_ov, right_ov, i):
    projs = []
    for k, wk in enumerate(penalty_states):
        Lk, Rk = left_ov[k].get(i - 1), right_ov[k].get(i + 2)
        pieces = [p for p in (Lk, wk.A(i), wk.A(i + 1), Rk) if p is not None]
        projs.append(contract_many(pieces))
    return projs


def _local_ground_state_penalized(L, Lbra, R, Rbra, H, ket, i, niter, projs, weight):
    """Same as _local_ground_state, but the matvec also adds
    weight*sum_k|proj_k><proj_k| (see this section's module-level note)."""
    matvec_H, order_in, shape, x0 = two_site_heff(L, Lbra, H, ket, i, R, Rbra)
    proj_flats = [p.transpose_to(order_in).reshape(-1) for p in projs]

    def matvec(v):
        out = matvec_H(v)
        for pf in proj_flats:
            # overlap stays on the device: it is multiplied straight back
            # into a device vector
            out = out + weight * bk.xp().vdot(pf, v) * pf
        return out

    dim = x0.size
    if dim <= 3:
        # dim <= 3, so this whole matrix is at most 3x3: build it with the
        # active namespace (matvec returns device arrays for a device
        # backend) and diagonalize it there rather than transferring.
        _xp = bk.xp()
        basis = bk.eye(dim)
        Hmat = _xp.column_stack([matvec(basis[:, k]) for k in range(dim)])
        w, v = _xp.linalg.eigh((Hmat + Hmat.conj().T) / 2)
        eval0, evec0 = w[0], v[:, 0]
    else:
        # The Sweeps schedule's own niter (e.g. ITensor's usual default of
        # 2) is deliberately small -- it relies on DMRG's many sweeps to
        # refine the state rather than fully converging each local
        # diagonalization -- but that only works once the sweep is already
        # close to converged; the OLD eigsh call floored this at
        # max(niter,200) (via maxiter=) for exactly this reason, and
        # dropping that floor when eigsh was replaced was a real
        # regression, confirmed directly: without it, larger local Hilbert
        # spaces (e.g. spin-3/2) and the excited-states penalty method
        # stopped converging to the true ground/target state (up to 4.7e-4
        # energy error on a case that should match ED to ~1e-14, and a
        # first-excited-state gap off by 2.7x). _lanczos_ground_state's own
        # early-stop-on-stable-Ritz-value logic means this floor costs
        # little in practice -- most bonds still exit in far fewer than
        # 200 iterations once bond dimension saturates.
        eval0, evec0 = _lanczos_ground_state(matvec, x0, niter=max(niter, 200))

    theta = ITensor(tuple(order_in), evec0.reshape(shape))
    return float(eval0.real), theta


def dmrg_excited(psi, H, penalty_states, weight, sweeps, quiet=True):
    """Two-site DMRG penalized against overlap with `penalty_states` (each
    weighted by `weight`) -- mirrors chain_session.h's
    dmrg(psi1,H_,wfs,sweeps,args-with-Weight) overload, used by
    excited_states() to find higher eigenstates one at a time. Mutates psi
    in place; returns the (penalized) objective's final local eigenvalue,
    which chain_session.h doesn't actually use as the reported energy
    either -- see excited_states() in chain.py, which recomputes
    <psi|H|psi> directly once converged.

    Convergence margin is genuinely thin here, more so than plain
    ground-state dmrg(): with no noise-term perturbation (see this
    module's docstring) to escape local minima of the penalized objective,
    a sweep over many random starting seeds at scale_lagrange=1.0 and a
    middling sweep count found a substantial fraction landing on the wrong
    (too-high) eigenvalue rather than slowly converging to the right one --
    not just occasional slow convergence, but real stationary points of
    the wrong energy. Doubling scale_lagrange and adding a few sweeps
    fixed every case checked *in that sweep*, matching standard DMRG
    practice that the penalty weight needs real margin above the
    bandwidth -- but this is not a universal fix. Confirmed directly on a
    different case (an 8-site transverse-field Ising chain, see
    examples/v2_VS_v3_excited_states_gap): the first-excited gap settled
    ~6e-4 off from the exact value, and this residual was *completely*
    insensitive to maxdim, sweep count (tested to 100, vs. the usual ~15-
    25), and scale_lagrange (tested 1x-100x bandwidth) -- a genuine
    stationary point the search cannot escape by throwing more of any of
    those knobs at it, for this Hamiltonian. So: more sweeps/weight is
    worth trying first (often sufficient, per the case above), but callers
    should treat a returned `fluctuations` entry that isn't small as a
    sign the search didn't actually converge and may not be fixable by
    parameter tuning alone -- only a real noise-term-style escape
    mechanism (not implemented here) would close this class of gap in
    general. This sensitivity is inherent to the algorithm as implemented,
    not a floating-point precision issue (verified directly against the
    JAX-vs-NumPy kernel comparison in kernels.py: a single matvec call
    agrees to ~1e-15 relative between the two, so tiny, ordinary floating-
    point differences between them are what's enough to occasionally tip
    an already-marginal case, not the root cause of a 6e-4 stationary
    point)."""
    n = psi.length()
    k_states = len(penalty_states)
    energy = None
    for sweep_i in range(sweeps.nsweep):
        maxdim, cutoff, noise, niter = sweeps.at(sweep_i)

        right_env = _all_right_environments(H, psi)
        right_ov = [_all_overlap_right(wk, psi) for wk in penalty_states]
        left_env = {0: (None, None)}
        left_ov = [{0: None} for _ in range(k_states)]
        for i in range(1, n):
            L, Lbra = left_env[i - 1]
            R, Rbra = right_env[i + 2]
            projs = _bond_projections(penalty_states, left_ov, right_ov, i)
            energy, theta = _local_ground_state_penalized(
                L, Lbra, R, Rbra, H, psi, i, niter, projs, weight)
            _apply_local_update(psi, i, theta, cutoff, maxdim, "right")
            left_env[i] = _extend_left(L, Lbra, H, psi, i)
            for k, wk in enumerate(penalty_states):
                left_ov[k][i] = _extend_overlap_left(left_ov[k][i - 1], wk, psi, i)

        right_env = {n + 1: (None, None)}
        right_ov = [{n + 1: None} for _ in range(k_states)]
        for i in range(n - 1, 0, -1):
            L, Lbra = left_env[i - 1]
            R, Rbra = right_env[i + 2]
            projs = _bond_projections(penalty_states, left_ov, right_ov, i)
            energy, theta = _local_ground_state_penalized(
                L, Lbra, R, Rbra, H, psi, i, niter, projs, weight)
            _apply_local_update(psi, i, theta, cutoff, maxdim, "left")
            right_env[i + 1] = _extend_right(R, Rbra, H, psi, i + 1)
            for k, wk in enumerate(penalty_states):
                right_ov[k][i + 1] = _extend_overlap_right(right_ov[k][i + 2], wk, psi, i + 1)
        left_env = {0: (None, None)}
        left_ov = [{0: None} for _ in range(k_states)]

        if not quiet:
            print("sweep {}: penalized energy = {}".format(sweep_i, energy))

    return energy
