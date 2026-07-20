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
from scipy.linalg import expm as _dense_expm

from .dmrg import (_all_left_environments, _all_right_environments,
                    _extend_left, _extend_right, one_site_heff, two_site_heff,
                    zero_site_heff)
from .mpsalgebra import _link_at
from .svd import svd
from .tensor import ITensor, commonIndex


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
    SVD truncation, since one-site TDVP must conserve bond dimension
    exactly), split it *losslessly* (QR-equivalent: svd with cutoff=0,
    maxdim=None) into a left-orthogonal Q kept at site i and a bond
    tensor C, then backward-evolve C (zero_site_heff, dmrg.py) and absorb
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
        Q, S, V, _spec = svd(A_new, ([left_link] if left_link else []) + [s_i],
                              cutoff=0.0, maxdim=None)
        psi.set_A(i, Q)
        new_link = commonIndex(Q, S)
        C = S * V
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
        U, S, V, _spec = svd(A_new, left_of_bond, cutoff=0.0, maxdim=None)
        psi.set_A(i, V)
        new_link = commonIndex(S, V)
        C = U * S
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
