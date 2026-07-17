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

    def matvec(v):
        t = ITensor(tuple(order_in), v.reshape(shape))
        for p in pieces:
            t = t * p
        return t.transpose_to(order_out).reshape(-1)

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

    def matvec(v):
        t = ITensor(tuple(order_in), v.reshape(shape))
        for p in pieces:
            t = t * p
        return t.transpose_to(order_out).reshape(-1)

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
