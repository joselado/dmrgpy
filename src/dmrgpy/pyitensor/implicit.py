"""Implicit (A-stable) imaginary-time steps of an MPS by two-site linear
solves.

Reference: J. P. Zima, E. M. Stoudenmire, S. R. White, O. Parcollet and
J. Kaye, "Fast Tensor Network Imaginary Time Evolution by Implicit
Stepping on Logarithmic Grids", arXiv:2606.02930.

One factor of a step maps psi to

    x = (1 + beta (H - s))^-1 (1 - beta (H - s)) psi,

with s a shift of the spectrum and beta a complex number. beta = dtau/2 is
the implicit trapezoid rule; the s-stage Gauss-Legendre rule is a product
of s such factors, beta = -dtau/q_j with q_j the roots of the Pade
denominator of exp(-z) (imaginarytime.py computes them). For Re(beta) > 0
and H - s >= 0 every factor has modulus at most one on every eigenstate,
whatever dtau is: this is the A-stability that lets the steps grow
exponentially on a logarithmic grid.

The linear system is solved the way ground-state DMRG solves its
eigenproblem: x is an MPS, swept two sites at a time, and at bond (i,i+1)
the two-site tensor theta of x solves the projected problem

    theta + beta (H_eff theta - s theta) = P[(1 + beta s) psi - beta H psi],

with H_eff the two-site effective Hamiltonian of x and P the projection
onto the tangent space of x at that bond (Galerkin, as in the paper's
Appendix B). Two choices differ from the paper, both to avoid a
compression the paper makes:

- No MPO for 1 + beta (H - s) is built: with x in mixed canonical gauge the
  projected identity is the identity, so the local operator is
  theta -> (1 - beta s) theta + beta H_eff theta, from the environments of
  H alone.
- The right-hand side is not formed as an MPS. P[psi] and P[H psi] come
  from the overlap environments <x|psi> and the mixed environments
  <x|H|psi>, so the apply-and-sum compression of b = (1 - beta(H - s)) psi
  is never made. The start of x is psi itself (fresh links), which is
  close to the solution once the state is near the ground state, where
  the steps are long.

The local system is not Hermitian for a complex beta, so it is solved by
restarted GMRES on the flat two-site vector; nothing truncates between
iterations, so the recurrence is exact. Every array goes through the
backend namespace, and only the small Hessenberg least-squares problem is
solved on the host.
"""

import numpy as np

from . import backend as bk
from .tensor import ITensor, contract_many, dag
from .mpsalgebra import _fresh_link_copy, inner
from .dmrg import (two_site_heff, _apply_local_update, _extend_left,
                   _extend_right, _all_right_environments)


def _fresh_copy(chain):
    """The same MPS with every Link index freshly minted, so that overlap
    environments between it and the original do not contract the links
    (ddmrg.py's _fresh_copy)."""
    cls = type(chain)
    out = cls(list(_fresh_link_copy(chain)))
    out.center = chain.center
    return out


def _overlap_left(L, x, b, i):
    piece = dag(x.A(i)) * b.A(i)
    return piece if L is None else L * piece


def _overlap_right(R, x, b, i):
    piece = dag(x.A(i)) * b.A(i)
    return piece if R is None else piece * R


def _local_overlap(Lov, Rov, b, i, order_in):
    """P[b] at bond (i,i+1), flattened in x's own (left link, s_i, s_j,
    right link) order."""
    pieces = [p for p in (Lov, b.A(i), b.A(i + 1), Rov) if p is not None]
    return contract_many(pieces).transpose_to(order_in).reshape(-1)


def _local_hpsi(L, Lbra, R, Rbra, H, b, i, order_out):
    """P[H b] at bond (i,i+1): the mixed environments <x|H|b> closed on b's
    two site tensors and H's two MPO tensors, in the bra order of the
    environments (x's links, primed site legs), which matches x's own
    order position by position."""
    pieces = [p for p in (L, H.A(i), H.A(i + 1), b.A(i), b.A(i + 1), R)
              if p is not None]
    return contract_many(pieces).transpose_to(order_out).reshape(-1)


def gmres(matvec, rhs, x0, tol=1e-10, krylov=10, maxiter=100):
    """Restarted GMRES for matvec(x) = rhs on flat vectors of the backend
    namespace. Stops when ||rhs - A x|| <= tol ||rhs|| or after maxiter
    matrix-vector products. Returns (x, relative residual)."""
    xp = bk.xp()
    bnorm = float(np.sqrt(abs(bk.scalar(xp.vdot(rhs, rhs)))))
    if bnorm == 0.0:
        return xp.zeros_like(rhs), 0.0
    x = x0
    r = rhs - matvec(x)
    beta = float(np.sqrt(abs(bk.scalar(xp.vdot(r, r)))))
    nmv = 1
    while beta > tol * bnorm and nmv < maxiter:
        m = min(krylov, maxiter - nmv)
        V = [r / beta]
        Hh = np.zeros((m + 1, m), dtype=complex)
        k_used = 0
        for k in range(m):
            w = matvec(V[k])
            nmv += 1
            for jj in range(k + 1):
                h = bk.scalar(xp.vdot(V[jj], w))
                Hh[jj, k] = h
                w = w - h * V[jj]
            hn = float(np.sqrt(abs(bk.scalar(xp.vdot(w, w)))))
            Hh[k + 1, k] = hn
            k_used = k + 1
            e1 = np.zeros(k_used + 1, dtype=complex)
            e1[0] = beta
            y, *_ = np.linalg.lstsq(Hh[:k_used + 1, :k_used], e1, rcond=None)
            res = float(np.linalg.norm(Hh[:k_used + 1, :k_used] @ y - e1))
            if hn == 0.0 or res <= tol * bnorm:
                break
            V.append(w / hn)
        upd = V[0] * y[0]
        for jj in range(1, k_used):
            upd = upd + V[jj] * y[jj]
        x = x + upd
        r = rhs - matvec(x)
        nmv += 1
        beta = float(np.sqrt(abs(bk.scalar(xp.vdot(r, r)))))
    return x, beta / bnorm


def implicit_factor(psi, H, beta, shift, cutoff, maxdim, nsweeps=1,
                    tol=1e-10, krylov=10, maxiter=100, x0=None):
    """x = (1 + beta(H - shift))^-1 (1 - beta(H - shift)) psi by two-site
    sweeps (see the module docstring). One sweep is a right and a left
    half sweep. The result has its orthogonality centre on site 1 and is
    not normalized. x0 optionally replaces psi as the start."""
    beta = complex(beta)
    a0 = 1.0 - beta * shift      # theta coefficient of the local operator
    b0 = 1.0 + beta * shift      # psi coefficient of the right-hand side
    x = _fresh_copy(x0 if x0 is not None else psi)
    n = x.length()
    x.position(1)

    def local(i, L, Lbra, R, Rbra, Lm, Lmb, Rm, Rmb, Lov, Rov):
        matvec_h, order_in, shape, xflat = two_site_heff(L, Lbra, H, x, i,
                                                          R, Rbra)
        s_i, s_j = [ind for ind in order_in if ind.hastags("Site")]
        order_out = (([Lmb] if Lmb is not None else []) +
                     [s_i.prime(1), s_j.prime(1)] +
                     ([Rmb] if Rmb is not None else []))
        rhs = (b0 * _local_overlap(Lov, Rov, psi, i, order_in)
               - beta * _local_hpsi(Lm, Lmb, Rm, Rmb, H, psi, i, order_out))

        def matvec(v):
            return a0 * v + beta * matvec_h(v)
        theta, _ = gmres(matvec, rhs, xflat, tol=tol, krylov=krylov,
                         maxiter=maxiter)
        return ITensor(tuple(order_in), theta.reshape(shape))

    for _ in range(nsweeps):
        right_env = _all_right_environments(H, x)
        right_mix = _all_right_environments(H, psi, bra=x)
        right_ov = {n + 1: None}
        for i in range(n, 1, -1):
            right_ov[i] = _overlap_right(right_ov[i + 1], x, psi, i)
        left_env = {0: (None, None)}
        left_mix = {0: (None, None)}
        left_ov = {0: None}
        for i in range(1, n):
            theta = local(i, *left_env[i - 1], *right_env[i + 2],
                          *left_mix[i - 1], *right_mix[i + 2],
                          left_ov[i - 1], right_ov[i + 2])
            _apply_local_update(x, i, theta, cutoff, maxdim, "right")
            left_env[i] = _extend_left(*left_env[i - 1], H, x, i)
            left_mix[i] = _extend_left(*left_mix[i - 1], H, psi, i, bra=x)
            left_ov[i] = _overlap_left(left_ov[i - 1], x, psi, i)
        right_env = {n + 1: (None, None)}
        right_mix = {n + 1: (None, None)}
        right_ov = {n + 1: None}
        for i in range(n - 1, 0, -1):
            theta = local(i, *left_env[i - 1], *right_env[i + 2],
                          *left_mix[i - 1], *right_mix[i + 2],
                          left_ov[i - 1], right_ov[i + 2])
            _apply_local_update(x, i, theta, cutoff, maxdim, "left")
            right_env[i + 1] = _extend_right(*right_env[i + 2], H, x, i + 1)
            right_mix[i + 1] = _extend_right(*right_mix[i + 2], H, psi,
                                             i + 1, bra=x)
            right_ov[i + 1] = _overlap_right(right_ov[i + 2], x, psi, i + 1)
    return x


def implicit_step(psi, H, betas, shift, cutoff, maxdim, nsweeps=1,
                  tol=1e-10, krylov=10, maxiter=100, shift_energy=False):
    """One imaginary-time step: the factors of `betas` applied in turn
    (implicit_factor), then the state normalized. shift_energy=True
    replaces `shift` by the energy <psi|H|psi>/<psi|psi> of the input
    state, the choice the paper makes for the Anderson model; a fixed
    shift at or below the ground-state energy keeps every factor
    A-stable."""
    if shift_energy:
        shift = (inner(psi, H, psi) / inner(psi, psi)).real
    x = psi
    for beta in betas:
        x = implicit_factor(x, H, beta, shift, cutoff, maxdim,
                            nsweeps=nsweeps, tol=tol, krylov=krylov,
                            maxiter=maxiter)
    x.normalize()
    return x
