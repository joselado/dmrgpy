"""Dynamical DMRG: the correction vector as a *variational* minimum
(Jeckelmann, Phys. Rev. B 66, 045114 (2002), cond-mat/0203500), instead of
a globally conjugate-gradient-solved linear system.

What this replaces, and why
---------------------------
`src/dmrgpy/cvm.py` solves

    M |x> = |b>,    M = (H - E0 - omega)^2 + eta^2,   |b> = -eta B|GS>

with a conjugate gradient built out of whole-MPS primitives: every CG
iterate is formed by applying an MPO to an MPS and adding MPS together,
and *each* of those steps compresses the result back down to `cvm_maxm`.
The compression is what breaks it. CG's residual recurrence assumes exact
arithmetic; a lossy projection after every operation puts a floor under
the reachable residual, and past that floor the recurrence does not merely
stall, it diverges -- which is why cvm.py carries `cvm_patience` /
`cvm_blowup` guards and returns the best iterate it happened to visit
rather than the one the loop ended on. Measured on a 16-site Heisenberg
chain at `cvm_maxm=40`, that floor sits at a residual of ~3e-2 against a
requested tolerance of 1e-5.

The variational formulation removes the floor at its source. Jeckelmann's
functional

    W(x) = <x|M|x> - 2 Re <b|x>

has its minimum exactly at M|x> = |b>, with W(x_min) = -<b|x_min>, and
minimizing it *inside* the MPS manifold -- one two-site tensor at a time,
sweeping, exactly as ground-state DMRG minimizes <psi|H|psi> -- means the
truncation is part of the variational ansatz rather than an error injected
after the fact. Each local problem is solved exactly (a small linear solve
on the two-site tensor), so there is no global recurrence to destabilize,
and the sequence of sweeps is monotonically non-increasing in W by
construction.

The value is second-order accurate on top of that: W is stationary at the
minimum, so an O(eps) error in x costs only O(eps^2) in W. Note that
cvm.py's own estimator already inherits this when A = B^dagger -- for that
case -Im<GS|A|x>/pi and -W_min/(pi*eta) are the same number, which is
Jeckelmann's point -- so the accuracy gained here comes from the better
x, not from a different formula applied to the same x.

Scope
-----
`itensor_version="python"` only. The two-site effective-operator
machinery this needs (environments, `two_site_heff`, the SVD-splitting
local update) already exists in `dmrg.py` and is reused verbatim; the
compiled backends would need the same solver written again in C++ against
ITensor's own environment classes, which is not done. `cvm.py` therefore
keeps its global CG as the default and for every other backend --
`Many_Body_Chain.cvm_solver` picks between them.

Off-diagonal correlators
------------------------
The variational principle is a statement about a real, bounded-below
quadratic form, which requires |b> to be the same operator on both sides:
it computes <GS|B^dagger (z-H)^-1 B|GS>. `cvm.py` supports a general
<GS|A (z-H)^-1 B|GS>, so this solver is used only when A = B^dagger (the
common case: any `name="XX"`-style diagonal correlator, and every i == j
correlator of a Hermitian operator). `cvm.py` checks and falls back.
"""

from . import backend as bk
from .tensor import ITensor, dag, contract_many
from .mpsalgebra import _fresh_link_copy, inner
from .dmrg import (two_site_heff, _apply_local_update, _extend_left,
                   _extend_right, _all_right_environments)


def _fresh_copy(chain):
    """The same MPS with every Link index freshly minted.

    Needed because this solver builds overlap environments between the
    ansatz `x` and the right-hand side `b`, and `x` is *initialized from*
    `b`: `dag(x.A(i)) * b.A(i)` on two chains that share Link Index
    objects would auto-contract the link legs as well as the physical one,
    silently collapsing the environment (the same trap `_fresh_link_copy`
    documents for inner()'s bra side)."""
    cls = type(chain)
    out = cls(list(_fresh_link_copy(chain)))
    out.center = chain.center
    return out


def _overlap_left(L, x, b, i):
    """One more site of the left environment <x_{<=i}|b_{<=i}>: free
    indices are (x's right link, b's right link)."""
    piece = dag(x.A(i)) * b.A(i)
    return piece if L is None else L * piece


def _overlap_right(R, x, b, i):
    piece = dag(x.A(i)) * b.A(i)
    return piece if R is None else piece * R


def _all_overlap_right(x, b):
    n = x.length()
    env = {n + 1: None}
    for i in range(n, 1, -1):
        env[i] = _overlap_right(env[i + 1], x, b, i)
    return env


def _local_rhs(Lov, Rov, b, i, order_in, shape):
    """The right-hand side |b> projected onto the two-site tangent space
    at bond (i,i+1), flattened in the same axis order `two_site_heff`
    hands back for the ansatz's own two-site tensor.

    Index bookkeeping: `Lov`/`Rov` carry one leg from x and one from b, so
    contracting them against b's two site tensors closes every b link and
    leaves exactly x's own (left link, s_i, s_j, right link) -- the space
    the local solve lives in."""
    pieces = [p for p in (Lov, b.A(i), b.A(i + 1), Rov) if p is not None]
    v = contract_many(pieces)
    return v.transpose_to(order_in).reshape(-1)


def _local_solve(matvec, rhs, x0, tol, maxiter):
    """Conjugate gradient on the *small* two-site system M_eff theta = rhs.

    CG is safe here in a way it is not for cvm.py's global loop: M_eff is
    Hermitian positive definite (it is the projection of
    (H-E0-omega)^2 + eta^2 onto a subspace) and, crucially, nothing
    truncates between iterations -- these are dense flat arrays, so the
    residual recurrence holds exactly and the usual convergence guarantee
    applies. Runs on whatever namespace the backend is using, so on a
    device the whole local solve stays on the device."""
    xp = bk.xp()
    x = x0
    r = rhs - matvec(x)
    p = r
    rs = bk.scalar(xp.vdot(r, r)).real
    b_norm = bk.scalar(xp.vdot(rhs, rhs)).real ** 0.5
    if b_norm == 0.0:
        return xp.zeros_like(rhs)
    for _ in range(maxiter):
        if rs ** 0.5 <= tol * b_norm:
            break
        Ap = matvec(p)
        denom = bk.scalar(xp.vdot(p, Ap)).real
        if denom == 0.0:
            break
        alpha = rs / denom
        x = x + alpha * p
        r = r - alpha * Ap
        rs_new = bk.scalar(xp.vdot(r, r)).real
        p = r + (rs_new / rs) * p
        rs = rs_new
    return x


def correction_vector(Mmpo, b, sweeps_maxdim, nsweeps, cutoff,
                      local_tol=1e-10, local_maxiter=200, x0=None,
                      quiet=True):
    """Minimize W(x) = <x|M|x> - 2 Re<b|x> over MPS of bond dimension at
    most `sweeps_maxdim`, by two-site sweeping.

    Returns (x, W_min). At the minimum M x = b, so W_min = -<b|x>; both
    are returned because the caller wants the *value* (second-order
    accurate) and, for the sum rule and for reporting, the vector.

    `x0` optionally seeds the ansatz -- unlike cvm.py's global CG, where
    warm-starting from a neighbouring frequency was measured to be
    actively harmful, a variational sweep can only ever be helped by a
    better starting point: every sweep is non-increasing in W regardless
    of where it starts.
    """
    x = _fresh_copy(x0 if x0 is not None else b)
    n = x.length()
    # Start the ansatz right-canonical so the first left-to-right sweep's
    # environments are the ones the local solve assumes.
    x.position(1)

    W = None
    for _ in range(nsweeps):
        right_env = _all_right_environments(Mmpo, x)
        right_ov = _all_overlap_right(x, b)
        left_env = {0: (None, None)}
        left_ov = {0: None}

        for i in range(1, n):
            L, Lbra = left_env[i - 1]
            R, Rbra = right_env[i + 2]
            matvec, order_in, shape, xflat = two_site_heff(
                L, Lbra, Mmpo, x, i, R, Rbra)
            rhs = _local_rhs(left_ov[i - 1], right_ov[i + 2], b, i,
                             order_in, shape)
            theta_flat = _local_solve(matvec, rhs, xflat,
                                      local_tol, local_maxiter)
            theta = ITensor(tuple(order_in), theta_flat.reshape(shape))
            _apply_local_update(x, i, theta, cutoff, sweeps_maxdim, "right")
            left_env[i] = _extend_left(L, Lbra, Mmpo, x, i)
            left_ov[i] = _overlap_left(left_ov[i - 1], x, b, i)

        right_env = {n + 1: (None, None)}
        right_ov = {n + 1: None}
        for i in range(n - 1, 0, -1):
            L, Lbra = left_env[i - 1]
            R, Rbra = right_env[i + 2]
            matvec, order_in, shape, xflat = two_site_heff(
                L, Lbra, Mmpo, x, i, R, Rbra)
            rhs = _local_rhs(left_ov[i - 1], right_ov[i + 2], b, i,
                             order_in, shape)
            theta_flat = _local_solve(matvec, rhs, xflat,
                                      local_tol, local_maxiter)
            theta = ITensor(tuple(order_in), theta_flat.reshape(shape))
            _apply_local_update(x, i, theta, cutoff, sweeps_maxdim, "left")
            right_env[i + 1] = _extend_right(R, Rbra, Mmpo, x, i + 1)
            right_ov[i + 1] = _overlap_right(right_ov[i + 2], x, b, i + 1)

        W = _functional(Mmpo, x, b)
        if not quiet:
            print("ddmrg sweep: W = %.12g" % W)

    return x, W


def _functional(Mmpo, x, b):
    """W(x) = <x|M|x> - 2 Re<b|x>, evaluated globally.

    Used as the sweep-to-sweep convergence monitor and as the estimator
    the caller turns into a spectral weight. Computed from whole-MPS
    primitives rather than from the last local bond's value, so it is the
    honest functional of the *whole* returned state."""
    quad = inner(x, Mmpo, x).real
    lin = inner(b, x).real
    return quad - 2.0 * lin
