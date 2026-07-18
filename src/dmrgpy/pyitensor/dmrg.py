"""Two-site DMRG ground-state solver.

Standard two-site sweep: build left/right environments of <psi|H|psi>,
diagonalize the local two-site effective Hamiltonian at each bond via
Lanczos (scipy.sparse.linalg.eigsh, applied as a LinearOperator -- the
2-site block is never densified beyond what eigsh itself needs), SVD-split
the result back into two site tensors with truncation, move the
orthogonality center, and repeat sweeps per the given Sweeps schedule.

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
"""

import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh

from . import kernels
from .mpsalgebra import _link_at
from .svd import svd
from .tensor import ITensor, dag
from .tensor import prime as _t_prime


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


def _extend_left(L, left_bra, H, ket, i):
    T = ket.A(i)
    bra_piece, _, right_bra = _relabel_bra_local(T, ket, i, left_bra, None)
    piece = bra_piece * H.A(i) * T
    new_L = piece if L is None else L * piece
    return new_L, right_bra


def _extend_right(R, right_bra, H, ket, i):
    T = ket.A(i)
    bra_piece, left_bra, _ = _relabel_bra_local(T, ket, i, None, right_bra)
    piece = bra_piece * H.A(i) * T
    new_R = piece if R is None else piece * R
    return new_R, left_bra


def _all_right_environments(H, ket):
    """{i: (R_tensor_or_None, dangling_bra_link_or_None)} for i = N+1..2."""
    n = ket.length()
    env = {n + 1: (None, None)}
    for i in range(n, 1, -1):
        R_next, bra_next = env[i + 1]
        env[i] = _extend_right(R_next, bra_next, H, ket, i)
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


def _local_ground_state(L, Lbra, R, Rbra, H, ket, i, niter):
    """Diagonalize the 2-site effective Hamiltonian at bond (i,i+1) for its
    lowest eigenpair. Returns (energy, theta_ITensor)."""
    matvec, order_in, shape, x0 = two_site_heff(L, Lbra, H, ket, i, R, Rbra)

    dim = x0.size
    if dim <= 3:
        # eigsh requires k < N; for a tiny effective space just diagonalize directly.
        basis = np.eye(dim, dtype=complex)
        Hmat = np.column_stack([matvec(basis[:, k]) for k in range(dim)])
        w, v = np.linalg.eigh((Hmat + Hmat.conj().T) / 2)
        eval0, evec0 = w[0], v[:, 0]
    else:
        Heff = LinearOperator((dim, dim), matvec=matvec, dtype=complex)
        ncv = min(dim, max(2 * niter + 1, 20))
        w, v = eigsh(Heff, k=1, which="SA", v0=x0, maxiter=max(niter, 200), ncv=ncv)
        eval0, evec0 = w[0], v[:, 0]

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


def dmrg(psi, H, sweeps, quiet=True):
    """Two-site ground-state DMRG. Mutates psi in place and returns the
    final energy -- mirrors ITensor's own dmrg(psi,H,sweeps,args) signature
    (minus the Args bag, which this package uses explicit kwargs for
    throughout, and minus the excited-state penalty overload -- see
    dmrg_excited() for that)."""
    n = psi.length()
    energy = None
    for sweep_i in range(sweeps.nsweep):
        maxdim, cutoff, noise, niter = sweeps.at(sweep_i)

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
        left_env = {0: (None, None)}

        if not quiet:
            print("sweep {}: energy = {}".format(sweep_i, energy))

    return energy


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
        p = wk.A(i) * wk.A(i + 1)
        p = p if Lk is None else Lk * p
        p = p if Rk is None else p * Rk
        projs.append(p)
    return projs


def _local_ground_state_penalized(L, Lbra, R, Rbra, H, ket, i, niter, projs, weight):
    """Same as _local_ground_state, but the matvec also adds
    weight*sum_k|proj_k><proj_k| (see this section's module-level note)."""
    matvec_H, order_in, shape, x0 = two_site_heff(L, Lbra, H, ket, i, R, Rbra)
    proj_flats = [p.transpose_to(order_in).reshape(-1) for p in projs]

    def matvec(v):
        out = matvec_H(v)
        for pf in proj_flats:
            out = out + weight * np.vdot(pf, v) * pf
        return out

    dim = x0.size
    if dim <= 3:
        basis = np.eye(dim, dtype=complex)
        Hmat = np.column_stack([matvec(basis[:, k]) for k in range(dim)])
        w, v = np.linalg.eigh((Hmat + Hmat.conj().T) / 2)
        eval0, evec0 = w[0], v[:, 0]
    else:
        Heff = LinearOperator((dim, dim), matvec=matvec, dtype=complex)
        ncv = min(dim, max(2 * niter + 1, 20))
        w, v = eigsh(Heff, k=1, which="SA", v0=x0, maxiter=max(niter, 200), ncv=ncv)
        eval0, evec0 = w[0], v[:, 0]

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
    fixed every case checked, matching standard DMRG practice that the
    penalty weight needs real margin above the bandwidth, but callers
    should treat a returned `fluctuations` entry that isn't small as a
    sign the search didn't actually converge, not just trust the energy.
    This sensitivity is inherent to the algorithm as implemented (verified
    directly against the JAX-vs-NumPy kernel comparison in kernels.py: a
    single matvec call agrees to ~1e-15 relative between the two, so tiny,
    ordinary floating-point differences between them are what's enough to
    occasionally tip a marginal case, not a correctness bug in either
    kernel)."""
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
