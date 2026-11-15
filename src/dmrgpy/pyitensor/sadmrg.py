"""State-averaged (multi-target) two-site DMRG: the k lowest eigenstates of a
Hermitian H in one sweep sequence, without an overlap penalty.

The ansatz is one MPS whose orthogonality center carries one extra index,
the state index, of dimension k (the "bundled" MPS of Baker, Foley and
Senechal, arXiv:2109.08181; block2's and Renormalizer's state-averaged
MPS). Every other tensor is shared by the k states, so <psi_i|psi_j> is the
overlap of their center slices and the environments of <psi|H|psi> are the
same for all of them. At each bond of a two-site sweep:

1. the k carried center vectors, contracted with the neighboring site, are
   the start block of a block Davidson solve for the k lowest eigenpairs of
   the two-site effective Hamiltonian (_block_davidson), preconditioned with
   its diagonal as the single-vector _davidson_ground_state of dmrg.py is;
2. the k Ritz vectors are stacked into one ITensor with the state index and
   split by the ordinary svd(), the state index grouped on the side the
   center moves to. The singular values of that stack are the eigenvalues of
   sum_j w_j theta_j theta_j^dag, so svd()'s truncation, which normalizes by
   the total weight, is the truncation of the state-averaged reduced density
   matrix (trace 1, equal weights 1/k unless `weights` says otherwise), with
   the same meaning of cutoff and maxdim as in dmrg();
3. the new center of each state is the projection of its own (unweighted)
   Ritz vector on the kept basis, dag(U)*theta or theta*dag(V).

After the last sweep the center, with the state index, is on site 1. The k
center slices are orthonormalized there and H is diagonalized in their span
(a Rayleigh-Ritz step with the one-site effective Hamiltonian, which is
exact for states that share every other tensor), and each state is cut out
as its own MPS with fresh link indices.

Layout. The state index never sits on the MPS itself during the sweep: the
working MPS holds the shared isometries and, at the center, the slice of
state 0 as a placeholder, while the k centers travel beside it as one
ITensor (`block`). This keeps two_site_heff, one_site_heff, _heff_diagonal
and the environment builders of dmrg.py usable unchanged, since they read
the index structure of the center from the MPS. It is Renormalizer's
layout (a list of per-root center tensors beside one MPS).

References and what was taken from each:

* Renormalizer (github.com/shuaigroup/Renormalizer, Apache-2.0),
  renormalizer/mps/gs.py and mp.py (_update_mps): the averaged density
  matrix with equal weights (sum over roots divided by their number), the
  per-root center as the projection on the kept basis, missing roots filled
  with random vectors, a dense eigh below a block size. Mirrored in
  structure, not copied: the truncation here is svd() of the weighted stack,
  which gives the same basis.
* PySCF lib/linalg_helper.py davidson1 (Apache-2.0): the multi-root
  Davidson shape, a restart onto the current Ritz vectors once the subspace
  exceeds a cap, and corrections only for roots not yet converged, dropped
  when linearly dependent. Departures: the stop is the per-root residual
  relative to the operator's scale, the rule of dmrg.py's
  _davidson_ground_state (PySCF tests an energy change and sqrt(tol) on the
  residual), and the preconditioner of root j uses its own Ritz value
  theta_j (PySCF passes the lowest one to every root), so that for k=1 the
  local solve is the Davidson iteration of dmrg.py's
  _davidson_ground_state.
* Baker, Foley, Senechal, arXiv:2109.08181 (bundled MPS, block Lanczos):
  the state index on the orthogonality center, moved through an SVD with the
  index on the side the center moves to, singular values pooled over the
  bundle.
* block2 documentation (state-averaged DMRG, nroots/weights; the code is
  GPL-3.0/LGPL-3.0 and was not read): solving more roots than needed when a
  requested k could cut a degenerate level, here `nguard`.

GPU. Every step runs on the array namespace of backend.py. What comes back
to the host: the Davidson Rayleigh matrix (at most SA_MAX_BLOCKS*k square),
the k residual norms per step, the singular values inside svd(), and the
final k x k Gram and Hamiltonian matrices. Random start vectors and the
first-sweep noise are drawn on the host from this module's own generator
and moved once, O(dim) per drawn vector. No array is written in place.

The local solves of one bond apply H_eff to the whole block at once
(SA_BATCHED_MATVEC): the block is one more free leg of the same contraction
plan kernels.make_matvec builds, so the k matvecs become one chain of
larger matrix products. SA_BATCHED_MATVEC = False applies the single-vector
matvec column by column, for comparison.

LAST_STATS holds the counters of the last dmrg_state_averaged() call:
matvecs (vector applications of H_eff, a block of b counting b), local
solves, dense solves, Davidson steps and restarts, the largest discarded
weight per sweep, the Ritz values of the last bond per sweep and the final
Rayleigh-Ritz values.
"""

import time
import warnings

import numpy as np

from . import backend as bk
from . import kernels
from .dmrg import (two_site_heff, one_site_heff, _heff_diagonal,
                   _extend_left, _extend_right, _all_right_environments)
from .index import Index
from .mpscontainer import MPS, _link_at
from .mpsalgebra import _fresh_link_copy, inner
from .svd import svd
from .tensor import ITensor, dag


# Davidson steps per local solve: SA_DAVIDSON_NITER on the early sweeps,
# SA_DAVIDSON_NITER_FINAL on the last SA_FULL_SWEEPS sweeps (the shape of
# dmrg.py's _penalized_budget), floored by the schedule's own niter. A step
# applies H_eff to the block of roots not yet converged.
SA_DAVIDSON_NITER = 4
SA_DAVIDSON_NITER_FINAL = 8
SA_FULL_SWEEPS = 2
# Per-root residual stop, relative to the operator's scale (the larger of
# max|theta| and max|diag H_eff|), and the Ritz-value stop: the values of
# _davidson_ground_state's ground-state solve.
SA_DAVIDSON_RTOL = 1e-12
SA_DAVIDSON_VTOL = 1e-14
# The Davidson subspace is restarted onto the current Ritz vectors once it
# would exceed SA_MAX_BLOCKS*nroots vectors. With SA_DAVIDSON_NITER=4 it
# never restarts.
SA_MAX_BLOCKS = 6
# A two-site space of dimension <= max(3, SA_DENSE_ROOT_FACTOR*nroots) is
# diagonalized densely: the dim <= 3 rule of _local_ground_state for one
# root, and the small bonds at the chain ends for many.
SA_DENSE_ROOT_FACTOR = 3
SA_BATCHED_MATVEC = True
# A candidate vector is dropped when what is left of it after two
# orthogonalization passes is below this times its norm (dmrg.py's value).
_DROP = 1e-10
# Final Rayleigh-Ritz: a direction of the center slices whose Gram
# eigenvalue is below this times the largest is dropped as linearly
# dependent (the states merged under truncation).
_GRAM_DROP = 1e-12

LAST_STATS = {}


def _new_stats():
    return dict(matvecs=0, solves=0, dense_solves=0, davidson_steps=0,
                restarts=0, truncerr=[], sweep_ritz=[], final_ritz=None,
                dropped_directions=0, nroots=None, solve_seconds=None)


def _sa_budget(sweep_i, nsweep):
    """Davidson step budget of sweep sweep_i (0-based) of nsweep."""
    if sweep_i >= nsweep - SA_FULL_SWEEPS:
        return SA_DAVIDSON_NITER_FINAL
    return SA_DAVIDSON_NITER


class _BlockHeff:
    """The two-site effective Hamiltonian of two_site_heff() at bond
    (i,i+1), applied to a block of vectors, a (dim, b) backend array.

    The pieces and the output index order are those two_site_heff() builds
    (L, W_i, W_{i+1}, R and the bra legs Lbra, s_i', s_j', Rbra); the block
    is one more free Index on the input and the output, so the plan is the
    same contraction with b carried through every matrix product. One plan
    per block width b, built on first use."""

    def __init__(self, L, Lbra, H, i, R, Rbra, order_in, shape, matvec, stats):
        s_i, s_j = [ind for ind in order_in if ind.hastags("Site")]
        self._pieces = [p for p in (L, H.A(i), H.A(i + 1), R) if p is not None]
        self._order_in = list(order_in)
        self._order_out = (([Lbra] if Lbra is not None else [])
                           + [s_i.prime(1), s_j.prime(1)]
                           + ([Rbra] if Rbra is not None else []))
        self._shape = tuple(shape)
        self._matvec = matvec
        self._plans = {}
        self._stats = stats
        self.dim = int(np.prod(self._shape))

    def apply(self, Xm):
        xp = bk.xp()
        b = int(Xm.shape[1])
        self._stats["matvecs"] += b
        if b == 1 or not SA_BATCHED_MATVEC:
            return xp.column_stack([self._matvec(Xm[:, j]) for j in range(b)])
        mv = self._plans.get(b)
        if mv is None:
            blk = Index(b, tags="Block")
            mv = kernels.make_matvec(self._pieces, self._order_in + [blk],
                                     self._shape + (b,),
                                     self._order_out + [blk])
            self._plans[b] = mv
        # C order: entry [p, c] of Xm is entry (p..., c) of the tensor with
        # indices order_in + [blk], and the output comes back the same way
        return mv(Xm.reshape(-1)).reshape(self.dim, b)


def _orthogonalize(t, B):
    """t minus its projection on the columns of B (orthonormal, or None),
    two passes, normalized; None when what is left is below _DROP times
    the norm of t, i.e. t lies in the span of B."""
    xp = bk.xp()
    n0 = float(bk.to_host(xp.linalg.norm(t)))
    if not n0 > 0.0:
        return None
    if B is not None:
        for _ in range(2):
            t = t - B @ (B.conj().T @ t)
    nt = float(bk.to_host(xp.linalg.norm(t)))
    if not nt > _DROP * n0:
        return None
    return t / nt


def _random_vector(dim, rng):
    return bk.asarray(rng.standard_normal(dim) + 1j * rng.standard_normal(dim))


def _noisy(x, noise, rng):
    """x plus a random vector of relative weight noise (norm
    sqrt(noise)*||x||), dmrg.py's _noisy_start with this module's own
    generator."""
    xp = bk.xp()
    r = _random_vector(int(x.size), rng)
    nx = xp.linalg.norm(x)
    return x + (np.sqrt(noise) * nx / xp.linalg.norm(r)) * r


def _start_block(starts, nb, noise, rng):
    """nb orthonormal start vectors: the carried ones (with noise mixed in
    when noise > 0), orthonormalized in order, then random vectors for the
    roots that are missing (the first bond of a fresh start carries one
    vector, a bond whose space was smaller than k carries fewer) or that
    truncation made linearly dependent."""
    xp = bk.xp()
    dim = int(starts[0].size)
    out = []
    for x in starts:
        if len(out) == nb:
            break
        if noise > 0.0:
            x = _noisy(x, noise, rng)
        t = _orthogonalize(x, xp.column_stack(out) if out else None)
        if t is not None:
            out.append(t)
    tries = 0
    while len(out) < nb:
        tries += 1
        if tries > 10 * nb + 10:
            raise RuntimeError(
                "state-averaged DMRG: could not complete %d orthonormal start "
                "vectors in a space of dimension %d" % (nb, dim))
        t = _orthogonalize(_random_vector(dim, rng),
                           xp.column_stack(out) if out else None)
        if t is not None:
            out.append(t)
    return out


def _dense_lowest(op, nb):
    """The nb lowest eigenpairs of a small H_eff, built densely by applying
    it to the identity (on the backend, as _local_ground_state does)."""
    xp = bk.xp()
    Hm = op.apply(bk.eye(op.dim))
    w, v = xp.linalg.eigh((Hm + Hm.conj().T) / 2)
    return np.asarray(bk.to_host(w[:nb]), dtype=float), v[:, :nb]


def _block_davidson(op, V0, diag, nroots, niter, rtol, vtol, stats):
    """The nroots lowest eigenpairs of the Hermitian H_eff `op` by block
    Davidson with the diagonal (Jacobi) preconditioner, started from the
    orthonormal vectors V0.

    Each step: Rayleigh-Ritz in the basis V (the m x m matrix V^dag H V on
    the host, m <= SA_MAX_BLOCKS*nroots), Ritz pairs (theta_j, x_j) of the
    nroots lowest, residuals r_j = H x_j - theta_j x_j. Root j counts as
    converged once ||r_j|| <= rtol*scale, scale the larger of max|theta|
    and max|diag|, the stop of _davidson_ground_state. The solve stops
    when every root has converged, when the budget of niter steps is spent,
    or when no Ritz value moved by more than vtol*scale since the previous
    step. Otherwise every unconverged root adds the correction
    (diag - theta_j)^-1 r_j (denominators clamped at 1e-12*scale),
    orthogonalized twice against the basis and the corrections of this
    step; when it lies in the basis the residual itself is taken, and when
    that lies in the basis too the root adds nothing. H_eff is applied to
    the new vectors as one block. When the basis would exceed the cap it is
    first replaced by the current Ritz vectors (thick restart, as PySCF's
    davidson1). For nroots=1 and no restart this is
    _davidson_ground_state step for step.

    Returns (theta, X): the nroots Ritz values on the host and the Ritz
    vectors, orthonormal columns of a (dim, nroots) backend array."""
    xp = bk.xp()
    dscale = float(bk.to_host(xp.max(xp.abs(diag)))) if diag.size else 0.0
    max_space = max(SA_MAX_BLOCKS * nroots, nroots + 1)
    Vm = xp.column_stack(V0)
    AVm = op.apply(Vm)
    theta_prev = None
    theta, X = None, None
    for step in range(niter):
        stats["davidson_steps"] += 1
        M = bk.to_host(Vm.conj().T @ AVm)
        w, s = np.linalg.eigh((M + M.conj().T) / 2)
        theta = np.asarray(w[:nroots], dtype=float)
        Sk = bk.asarray(s[:, :nroots])
        X = Vm @ Sk
        AX = AVm @ Sk
        Rm = AX - X * bk.asarray(theta)[None, :]
        rn = np.asarray(bk.to_host(xp.linalg.norm(Rm, axis=0)), dtype=float)
        scale = max(float(np.max(np.abs(theta))), dscale)
        todo = [j for j in range(nroots) if rn[j] > rtol * scale]
        if step == niter - 1 or not todo:
            break
        if (vtol is not None and theta_prev is not None
                and float(np.max(np.abs(theta - theta_prev))) < vtol * scale):
            break
        theta_prev = theta
        if Vm.shape[1] + len(todo) > max_space:
            Vm, AVm = X, AX
            stats["restarts"] += 1
        floor = 1e-12 * scale if scale > 0 else 1.0
        new = []
        for j in todo:
            B = Vm if not new else xp.concatenate([Vm, xp.column_stack(new)], axis=1)
            den = diag - theta[j]
            den = xp.where(xp.abs(den) < floor, floor, den)
            t = _orthogonalize(Rm[:, j] / den, B)
            if t is None:
                # the preconditioned correction lies in the basis (H nearly
                # diagonal): take the residual, the Lanczos direction
                t = _orthogonalize(Rm[:, j], B)
            if t is not None:
                new.append(t)
        if not new:
            break
        T = xp.column_stack(new)
        Vm = xp.concatenate([Vm, T], axis=1)
        AVm = xp.concatenate([AVm, op.apply(T)], axis=1)
    X = X / xp.linalg.norm(X, axis=0)[None, :]
    return theta, X


def _unit(ind, j):
    """The unit vector e_j on Index ind, as an ITensor: contracting a block
    with it cuts out member j."""
    e = np.zeros(ind.dim)
    e[j] = 1.0
    return ITensor((ind,), e)


def _bond_update(work, H, i, L, Lbra, R, Rbra, block, direction, nroots,
                 sqrt_w, maxdim, cutoff, budget, noise, rng, stats):
    """One state-averaged local update at bond (i,i+1). `block` is (C, st),
    the carried centers with their state index st, on site i for a right
    move and on site i+1 for a left one, or None on the very first bond of a
    fresh start (the start is then work's own two-site tensor). Splits the
    result into work (isometry plus the state-0 placeholder) and returns
    (theta, (C_new, st_new), truncerr)."""
    matvec, order_in, shape, x0 = two_site_heff(L, Lbra, H, work, i, R, Rbra)
    order_in = list(order_in)
    shape = tuple(shape)
    dim = int(x0.size)
    nb = min(nroots, dim)
    if block is None:
        starts = [x0]
    else:
        C, st = block
        T = C * work.A(i + 1) if direction == "right" else work.A(i) * C
        Xs = T.transpose_to(order_in + [st]).reshape(dim, st.dim)
        starts = [Xs[:, j] for j in range(st.dim)]
    V0 = _start_block(starts, nb, noise, rng)

    op = _BlockHeff(L, Lbra, H, i, R, Rbra, order_in, shape, matvec, stats)
    stats["solves"] += 1
    if dim <= max(3, SA_DENSE_ROOT_FACTOR * nroots):
        stats["dense_solves"] += 1
        theta, X = _dense_lowest(op, nb)
    else:
        diag = _heff_diagonal(L, Lbra, H, i, R, Rbra, order_in)
        theta, X = _block_davidson(op, V0, diag, nb, budget, SA_DAVIDSON_RTOL,
                                   SA_DAVIDSON_VTOL, stats)

    st_new = Index(nb, tags="State")
    inds = tuple(order_in) + (st_new,)
    theta_u = ITensor(inds, X.reshape(shape + (nb,)))
    if sqrt_w is None:
        theta_w = theta_u
    else:
        Xw = X * bk.asarray(sqrt_w[:nb])[None, :]
        theta_w = ITensor(inds, Xw.reshape(shape + (nb,)))

    left_link = _link_at(work, i, i - 1)
    s_i = [ind for ind in order_in if ind.hastags("Site")][0]
    left_inds = ([left_link] if left_link is not None else []) + [s_i]
    if direction == "right":
        U, S, V, spec = svd(theta_w, left_inds, cutoff=cutoff, maxdim=maxdim)
        C_new = dag(U) * theta_u
        work.set_A(i, U)
        work.set_A(i + 1, C_new * _unit(st_new, 0))
        work.center = i + 1
    else:
        U, S, V, spec = svd(theta_w, left_inds + [st_new], cutoff=cutoff,
                            maxdim=maxdim)
        C_new = theta_u * dag(V)
        work.set_A(i + 1, V)
        work.set_A(i, C_new * _unit(st_new, 0))
        work.center = i
    return theta, (C_new, st_new), float(spec.truncerr)


def _extract_states(work, H, block, R2, k, stats):
    """The k states of the final center (site 1): orthonormalize the center
    slices and diagonalize H in their span with the one-site effective
    Hamiltonian at site 1 (exact, since every other tensor is shared and
    right-orthonormal), then cut out each state as its own MPS with fresh
    link indices. Returns the list of MPS, lowest Ritz value first."""
    xp = bk.xp()
    C, st = block
    R, Rbra = R2
    mv1, order1, shape1, _ = one_site_heff(None, None, H, work, 1, R, Rbra)
    order1 = list(order1)
    shape1 = tuple(shape1)
    nb = st.dim
    Cm = C.transpose_to(order1 + [st]).reshape(-1, nb)
    HC = xp.column_stack([mv1(Cm[:, j]) for j in range(nb)])
    stats["matvecs"] += nb
    G = bk.to_host(Cm.conj().T @ Cm)
    Hm = bk.to_host(Cm.conj().T @ HC)
    G = (G + G.conj().T) / 2
    Hm = (Hm + Hm.conj().T) / 2
    g, Q = np.linalg.eigh(G)
    keep = g > _GRAM_DROP * float(np.max(g))
    stats["dropped_directions"] = int(np.sum(~keep))
    Wm = Q[:, keep] / np.sqrt(g[keep])[None, :]
    e, Y = np.linalg.eigh(Wm.conj().T @ Hm @ Wm)
    Z = Wm @ Y
    Cn = Cm @ bk.asarray(Z)
    stats["final_ritz"] = [float(x) for x in e]
    nout = min(k, int(Cn.shape[1]))
    if nout < k:
        warnings.warn(
            "state-averaged DMRG: only %d of the %d requested states are "
            "linearly independent at the end of the sweeps; raise maxm"
            % (nout, k))
    states = []
    for j in range(nout):
        tmp = work.copy()
        tmp.set_A(1, ITensor(tuple(order1), Cn[:, j].reshape(shape1)))
        psi_j = MPS(_fresh_link_copy(tmp))
        psi_j.center = 1
        psi_j.normalize()
        states.append(psi_j)
    return states


def dmrg_state_averaged(psi, H, k, sweeps, quiet=True, nguard=0, weights=None,
                        noisy_start=True, seed=0):
    """The k lowest eigenstates of a Hermitian MPO H by state-averaged
    two-site DMRG (see this module's docstring), from the start MPS psi
    (not modified), with the Sweeps schedule's maxdim/cutoff per sweep.

    Returns (energies, states): the energies <psi_j|H|psi_j> of the k
    returned states, ascending, and the states as separate normalized MPS,
    mutually orthonormal.

    nguard solves and carries nguard more roots than are returned, all in
    the averaged density matrix: when k would cut a degenerate level (k=2
    on a singlet ground state with a triplet above), the k-th root is an
    arbitrary member of that level, a different one at each bond, and the
    shared basis follows it; the guard roots complete the level. weights,
    when given, are k+nguard nonnegative weights of the averaged density
    matrix (normalized here); None is equal weights. noisy_start mixes the
    first sweep's noise into every local start of that sweep, as dmrg()
    does on a fresh start. seed seeds this call's own generator (random
    fill vectors and noise), so the same call gives the same result.

    k=1 is ordinary two-site DMRG (one root, the Davidson iteration of
    _davidson_ground_state, the same SVD split), but not dmrg() bit for
    bit: the last SA_FULL_SWEEPS sweeps take SA_DAVIDSON_NITER_FINAL steps
    where dmrg() takes DAVIDSON_NITER, the first-sweep noise comes from
    this call's own generator, and the energy is measured on the returned
    state rather than read off the last local solve.

    LAST_STATS["solve_seconds"] is the wall time of the sweeps and the
    extraction, without the caller's own measurements."""
    global LAST_STATS
    k = int(k)
    nguard = int(nguard)
    if k < 1 or nguard < 0:
        raise ValueError("dmrg_state_averaged: need k >= 1 and nguard >= 0, "
                         "got k=%d, nguard=%d" % (k, nguard))
    n = psi.length()
    if n < 2:
        raise ValueError("dmrg_state_averaged: a two-site sweep needs at "
                         "least 2 sites, got %d" % n)
    nroots = k + nguard
    sqrt_w = None
    if weights is not None:
        w = np.asarray(weights, dtype=float).reshape(-1)
        if len(w) != nroots or np.any(w < 0) or not np.sum(w) > 0:
            raise ValueError("dmrg_state_averaged: weights must be %d "
                             "nonnegative numbers, not all zero" % nroots)
        sqrt_w = np.sqrt(w / np.sum(w))
    stats = _new_stats()
    LAST_STATS = stats
    rng = np.random.default_rng(seed)
    t_start = time.perf_counter()
    stats["nroots"] = nroots

    work = psi.copy()
    work.position(1)
    block = None
    right_env = _all_right_environments(H, work)
    for sweep_i in range(sweeps.nsweep):
        maxdim, cutoff, noise, niter = sweeps.at(sweep_i)
        budget = max(niter, _sa_budget(sweep_i, sweeps.nsweep))
        nz = noise if (noisy_start and sweep_i == 0) else 0.0
        terr = 0.0
        left_env = {0: (None, None)}
        for i in range(1, n):
            L, Lbra = left_env[i - 1]
            R, Rbra = right_env[i + 2]
            theta, block, te = _bond_update(work, H, i, L, Lbra, R, Rbra, block,
                                            "right", nroots, sqrt_w, maxdim,
                                            cutoff, budget, nz, rng, stats)
            terr = max(terr, te)
            left_env[i] = _extend_left(L, Lbra, H, work, i)
        right_env = {n + 1: (None, None)}
        for i in range(n - 1, 0, -1):
            L, Lbra = left_env[i - 1]
            R, Rbra = right_env[i + 2]
            theta, block, te = _bond_update(work, H, i, L, Lbra, R, Rbra, block,
                                            "left", nroots, sqrt_w, maxdim,
                                            cutoff, budget, nz, rng, stats)
            terr = max(terr, te)
            right_env[i + 1] = _extend_right(R, Rbra, H, work, i + 1)
        stats["truncerr"].append(terr)
        stats["sweep_ritz"].append([float(x) for x in theta])
        if not quiet:
            print("sweep {}: Ritz values = {}, discarded weight = {:.2e}".format(
                sweep_i, np.array2string(np.asarray(theta), precision=10), terr))

    states = _extract_states(work, H, block, right_env[2], k, stats)
    energies = [float(inner(s, H, s).real) for s in states]
    stats["solve_seconds"] = time.perf_counter() - t_start
    order = np.argsort(energies, kind="stable")
    return [energies[j] for j in order], [states[j] for j in order]
