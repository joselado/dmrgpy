"""Two-site TDVP real-time evolution.

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
"""

import numpy as np
from scipy.linalg import expm as _dense_expm

from .dmrg import (_all_left_environments, _all_right_environments,
                    _extend_left, _extend_right, one_site_heff, two_site_heff)
from .mpsalgebra import _link_at
from .svd import svd
from .tensor import ITensor


def _lanczos_expm_multiply(matvec, v0, coeff, niter=30, tol=1e-12):
    """expm(coeff * A) @ v0 for a Hermitian linear operator A (given as a
    matvec function) via a Krylov (Lanczos) subspace of dimension up to
    `niter`: build an orthonormal basis {v0, A v0, A^2 v0, ...} via the
    Lanczos recursion (with full reorthogonalization -- the subspace is
    small enough, niter ~ tens, that this costs nothing and buys real
    numerical stability), project A onto it as a small tridiagonal matrix
    T, exponentiate T exactly (cheap, dense), and map back. Standard
    Krylov-propagator technique for real-time quantum dynamics; this is
    what tdvp_step()'s "niter=50 bounds the Lanczos iterations" comment in
    mpscpp3/chain_session.h refers to on the compiled-backend side."""
    beta0 = np.linalg.norm(v0)
    if beta0 == 0:
        return v0.copy()
    q = v0 / beta0
    qs = [q]
    alphas = []
    betas = []
    w = matvec(q)
    alpha = np.vdot(q, w).real
    alphas.append(alpha)
    w = w - alpha * q
    m = min(niter, v0.size)
    for _ in range(1, m):
        beta = np.linalg.norm(w)
        if beta < tol:
            break
        betas.append(beta)
        q_new = w / beta
        qs.append(q_new)
        w = matvec(q_new)
        alpha = np.vdot(q_new, w).real
        alphas.append(alpha)
        w = w - alpha * q_new - beta * qs[-2]
        for qk in qs[:-1]:
            w = w - np.vdot(qk, w) * qk

    k = len(alphas)
    T = np.zeros((k, k), dtype=complex)
    for idx in range(k):
        T[idx, idx] = alphas[idx]
    for idx in range(k - 1):
        T[idx, idx + 1] = betas[idx]
        T[idx + 1, idx] = betas[idx]
    Q = np.column_stack(qs[:k])
    expT = _dense_expm(coeff * T)
    return beta0 * (Q @ expT[:, 0])


def _evolve_two_site(L, Lbra, H, ket, i, R, Rbra, tau, niter):
    matvec, order_in, shape, x0 = two_site_heff(L, Lbra, H, ket, i, R, Rbra)
    evolved = _lanczos_expm_multiply(matvec, x0, -1j * tau, niter=niter)
    return ITensor(tuple(order_in), evolved.reshape(shape))


def _evolve_one_site(L, Lbra, H, ket, i, R, Rbra, tau, niter):
    matvec, order_in, shape, x0 = one_site_heff(L, Lbra, H, ket, i, R, Rbra)
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


def tdvp_step(psi, H, dt, cutoff, maxdim, niter=50):
    """One real-time step exp(-i*dt*H) via two-site TDVP: a left-to-right
    half-sweep evolving by dt/2, then a right-to-left half-sweep evolving
    by another dt/2 -- mirrors mpscpp3/chain_session.h's tdvp_step()
    (NumCenter=2 there; this module doesn't implement one-site TDVP or the
    global subspace-expansion machinery, matching that file's own
    comment on why it's unnecessary here). Mutates psi in place."""
    tau = dt / 2.0
    _half_sweep_lr(psi, H, tau, cutoff, maxdim, niter)
    _half_sweep_rl(psi, H, tau, cutoff, maxdim, niter)
    return psi
