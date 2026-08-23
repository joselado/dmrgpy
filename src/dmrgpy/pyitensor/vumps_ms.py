"""Sequential multi-site VUMPS: a uniform MPS with a unit cell of ANY number
of sites, at a cost LINEAR rather than exponential in the cell size.

== Why this module exists ==

`vumps.py` supports a multi-site unit cell by *grouping*: `_group_automaton`
folds the cell's `n_uc` sites into one supersite of dimension
`d_g = prod(d_p)`, and the single-site algorithm then runs on that. This is
exact, and it is why that module rejects `n_uc > 2` -- a 4-site spinful cell
would be `d_g = 256`, and the two-site objects built on it are `d_g**2`.

That exponential blow-up is a known trap rather than an implementation
detail. Nietner, Vanhecke, Verstraete, Eisert and Vanderstraeten
(arXiv:2003.01142) put it directly: "the cost of a naive application of the
VUMPS algorithm would scale exponentially with the size of the unit cell. It
is therefore desirable to work with a variational algorithm that directly
takes into account the non-trivial structure of the unit cell", and their
stated key property is "a computational effort that scales linearly rather
than exponentially in the size of the unit cell". The multi-site VUMPS of
Zauner-Stauber, Vanderstraeten, Fishman, Verstraete and Haegeman (PRB 97,
045145 (2018), arXiv:1701.07035) is the algorithm that does this, and it is
what production codes implement (TeNPy's `SingleSiteVUMPSEngine` /
`TwoSiteVUMPSEngine` are `Sweep` subclasses that sweep over the MPS unit
cell; MPSKit likewise).

This module implements that: the state is a LIST of per-site tensors
`AL[n]`, `AR[n]`, `C[n]`, and one iteration sweeps over the cell solving a
one-site eigenproblem `H_AC[n]` at each site and a zero-site eigenproblem
`H_C[n]` at each bond. Nothing is ever grouped, so the cost per iteration is
`n_uc` times a single-site cost.

== What it generalizes, beyond the site count ==

`vumps.py`'s environments are specialized to a *reach-1* automaton: it keeps
one accumulated `GL` matrix, one `GR`, and a list of one-site-away bond
channels, and `idmrg_excitations._check_reach_one` rejects anything longer.
That specialization is exactly what grouping buys -- any coupling inside the
cell becomes on-site once the cell is one supersite.

Ungrouped, that is no longer true (a 3-site cell may couple sites 0 and 2),
so this module keeps the FULL channel-resolved environment instead:
`GL[n]` and `GR[n]` are `(Dw, D, D)` arrays, one `(D,D)` matrix per MPO
channel. That is the textbook formulation, it subsumes the reach-1 one, and
it removes the reach restriction rather than trading it for another.

The automaton convention is `idmrg._build_periodic_mpo`'s, unchanged:
`W[chan_l, s_in, s_out, chan_r]` with channel `_S_IDX` = "no term started
yet" and `_F_IDX` = "term finished / accumulating", so `W[S,:,:,S]` and
`W[F,:,:,F]` are identities, `W[S,:,:,F]` is the on-site term, and the
pending channels in between carry partially-applied terms. In that language
the environments are

    GL[S] = I                      (nothing applied to the left)
    GL[pending]                    (finite: the pending block is nilpotent)
    GL[F]                          (the accumulated Hamiltonian: geometric
                                    series, hence a regularized solve)

and mirrored on the right, where the roles of S and F swap (`GR[F] = I`).
Reading `vumps.py`'s own `_h_c_action` in this language shows it is the
`n_uc = 1`, reach-1 special case of the formulas here: its `_cap_right(C,GR)`
is the `S` channel, its `_cap_left(GL,C)` is `F`, and its bond terms are the
pending channels. That correspondence is what makes the grouped module a
usable oracle -- `tests/test_vumps_multisite.py` checks the two against each
other at `n_uc = 1` and `2`, where both are valid.

== Index conventions ==

Taken verbatim from the existing, validated helpers rather than re-derived
(`idmrg_excitations._op_transfer_matrix`, `idmrg._apply_transfer`,
`_apply_transfer_from_left`), because getting them subtly wrong is the
classic way to produce a plausible-looking wrong energy:

* site tensors are `(left, phys, right)`;
* an operator matrix is `M[in, out]`, applied to the ket only
  (`idmrg_excitations._apply_op_ket`);
* a left environment is `X[ket, bra]` and pushes forward as
  `X'[r,R] = sum X[l,L] (M.A)[l,p,r] conj(A)[L,p,R]`;
* a right environment is `Y[ket, bra]` and pushes backward as
  `Y'[l,L] = sum (M.A)[l,p,r] conj(A)[L,p,R] Y[r,R]`.
"""

import numpy as np
from scipy.linalg import polar as _polar

from . import idmrg
from . import idmrg_excitations as idmrg_exc
from .dmrg import _lanczos_ground_state
from .idmrg_excitations import _F_IDX, _S_IDX


# -- transfer primitives -----------------------------------------------------

def _push_left(X, A, M=None):
    """One site forward for a left environment: `X[l,L] -> X'[r,R]`, with
    operator `M` (M[in,out]) applied to the ket. Staged as two tensordots
    so no `(D,D,D,D)` transfer matrix is ever materialized -- `O(D^2 d D)`
    work instead of `O(D^4)`, the same avoidance `vumps.py`'s own
    `_precompute_bond_environments` documents."""
    B = idmrg_exc._apply_op_ket(M, A)
    # einsum('lL,lpr->Lpr', X, B) then einsum('Lpr,LpR->rR', ., conj(A))
    T = np.tensordot(X, B, axes=([0], [0]))
    return np.tensordot(T, np.conj(A), axes=([0, 1], [0, 1]))


def _push_right(Y, A, M=None):
    """One site backward for a right environment: `Y[r,R] -> Y'[l,L]`."""
    B = idmrg_exc._apply_op_ket(M, A)
    # einsum('lpr,rR->lpR', B, Y) then einsum('lpR,LpR->lL', ., conj(A))
    T = np.tensordot(B, Y, axes=([2], [0]))
    return np.tensordot(T, np.conj(A), axes=([1, 2], [1, 2]))


def _push_cell_left(G, A_list, W_list):
    """A whole unit cell forward for a channel-resolved left environment
    `G[a][l,L]`: for each site, `G'[b] = sum_a push(G[a], W[a,:,:,b])`.
    Only the channel pairs `(a,b)` whose `W` slice is actually nonzero are
    contracted -- the automaton is sparse (an upper-triangular block
    structure), so skipping the zeros is what keeps this linear."""
    for A, W in zip(A_list, W_list):
        Dw = W.shape[0]
        new = np.zeros((Dw,) + G.shape[1:], dtype=complex)
        for a in range(Dw):
            if not G[a].any():
                continue
            for b in range(Dw):
                slab = W[a, :, :, b]
                if not slab.any():
                    continue
                new[b] += _push_left(G[a], A, slab)
        G = new
    return G


def _push_cell_right(G, A_list, W_list):
    """Mirror of `_push_cell_left`: a whole cell backward, walking the
    sites in reverse, `G'[a] = sum_b push(G[b], W[a,:,:,b])`."""
    for A, W in zip(reversed(A_list), reversed(W_list)):
        Dw = W.shape[0]
        new = np.zeros((Dw,) + G.shape[1:], dtype=complex)
        for b in range(Dw):
            if not G[b].any():
                continue
            for a in range(Dw):
                slab = W[a, :, :, b]
                if not slab.any():
                    continue
                new[a] += _push_right(G[b], A, slab)
        G = new
    return G


# -- environments ------------------------------------------------------------

def _nilpotent_channels_left(A_list, W_list, D, Dw):
    """`(G, source)`: the channel-resolved left environment at the cell's
    left edge for every channel EXCEPT `F`, plus the `F` content one cell
    traversal produces from it.

    `G[S] = I` is fixed by definition, and the pending block is nilpotent
    (a Hamiltonian term has finite reach, so a partially-applied term must
    close within a bounded number of sites), which is why this converges by
    plain iteration instead of needing a solve: each sweep propagates the
    pending content one cell further, and after the nilpotency degree is
    exhausted nothing changes. `Dw + 1` sweeps is a hard bound on that
    degree; the loop exits as soon as it stops moving, which for the
    nearest-neighbour models this is usually run on is after one or two.

    `F` alone is left to `_solve_left_environment_cell`, because `F` is the
    one channel that feeds back into itself (`W[F,:,:,F] = I`) -- summing
    that is a geometric series over the whole half-infinite chain, not a
    finite propagation."""
    G = np.zeros((Dw, D, D), dtype=complex)
    G[_S_IDX] = np.eye(D, dtype=complex)
    for _ in range(Dw + 1):
        new = _push_cell_left(G, A_list, W_list)
        new[_S_IDX] = np.eye(D, dtype=complex)
        new[_F_IDX] = 0.0
        if np.allclose(new, G, atol=1e-14, rtol=0.0):
            G = new
            break
        G = new
    source = _push_cell_left(G, A_list, W_list)[_F_IDX]
    return G, source


def _nilpotent_channels_right(A_list, W_list, D, Dw):
    """Mirror of `_nilpotent_channels_left`, with `S` and `F` swapping
    roles: on the right it is `F` that means "nothing applied yet"
    (`G[F] = I`) and `S` that accumulates."""
    G = np.zeros((Dw, D, D), dtype=complex)
    G[_F_IDX] = np.eye(D, dtype=complex)
    for _ in range(Dw + 1):
        new = _push_cell_right(G, A_list, W_list)
        new[_F_IDX] = np.eye(D, dtype=complex)
        new[_S_IDX] = 0.0
        if np.allclose(new, G, atol=1e-14, rtol=0.0):
            G = new
            break
        G = new
    source = _push_cell_right(G, A_list, W_list)[_S_IDX]
    return G, source


def _cell_transfer_action_left(X, AL_list):
    """`X -> T_cell[X]`, the plain (operator-free) transfer of the whole
    cell, applied to a left environment."""
    for A in AL_list:
        X = _push_left(X, A)
    return X


def _cell_transfer_action_right(Y, AR_list):
    for A in reversed(AR_list):
        Y = _push_right(Y, A)
    return Y


def _solve_left_environment_cell(AL_list, r_AL, source, e, D):
    """`GL[F]`: solves `(I - T_cell + P)[G] = source - e*I` with
    `P[X] = I*trace(conj(r_AL) @ X)`.

    This is `vumps.py`'s `_solve_left_environment` with the single-site
    transfer replaced by the whole cell's, and it inherits that function's
    two hard-won requirements verbatim: the conjugate on `r_AL` is
    required, not optional (it is the fixed point of the *adjoint* map, and
    dropping it reproduces a plausible but wrong energy density), and `e`
    MUST be this side's own `e_L` rather than any average with the right
    side's -- a shared regularization target leaves the right-hand side
    with a nonzero component along the null direction `P` exists to remove,
    which showed up as an exact period-2 limit cycle. See that function's
    docstring for both."""
    I = np.eye(D, dtype=complex)

    def action(X):
        return (X - _cell_transfer_action_left(X, AL_list)
                + I * np.trace(r_AL.conj() @ X))

    return idmrg_exc._solve_linear_map(D, action, source - e * I)


def _solve_right_environment_cell(AR_list, l_AR, source, e, D):
    """Mirror of `_solve_left_environment_cell`."""
    I = np.eye(D, dtype=complex)

    def action(X):
        return (X - _cell_transfer_action_right(X, AR_list)
                + I * np.trace(l_AR.conj() @ X))

    return idmrg_exc._solve_linear_map(D, action, source - e * I)


def _cell_fixed_points(AL_list, AR_list, D):
    """`(r_AL, l_AR)`: the dominant right fixed point of the whole cell's
    AL-transfer and the dominant left fixed point of its AR-transfer, both
    Hermitized. Obtained matrix-free by powering the cell action rather
    than by building the `(D^2, D^2)` cell transfer matrix, so this stays
    linear in the cell length too."""
    from scipy.sparse.linalg import LinearOperator, eigs

    def _dominant(action, D):
        n = D * D
        if n <= 4:
            mat = np.zeros((n, n), dtype=complex)
            basis = np.eye(n, dtype=complex)
            for k in range(n):
                mat[:, k] = action(basis[:, k].reshape(D, D)).reshape(-1)
            w, v = np.linalg.eig(mat)
            k = int(np.argmax(np.abs(w)))
            return v[:, k].reshape(D, D)
        op = LinearOperator((n, n), dtype=complex,
                             matvec=lambda x: action(x.reshape(D, D)).reshape(-1))
        w, v = eigs(op, k=1, which='LM', maxiter=5000, tol=1e-13)
        return v[:, 0].reshape(D, D)

    r_AL = _dominant(lambda X: _cell_transfer_action_right(X, AL_list), D)
    l_AR = _dominant(lambda X: _cell_transfer_action_left(X, AR_list), D)
    r_AL = (r_AL + r_AL.conj().T) / 2
    l_AR = (l_AR + l_AR.conj().T) / 2
    # Fix the scale the eigensolver leaves arbitrary: these close a
    # normalized state, so their trace against the other side's exact
    # (identity) fixed point must be 1.
    tr_r = np.trace(r_AL)
    if abs(tr_r) > 1e-300:
        r_AL = r_AL / tr_r
    tr_l = np.trace(l_AR)
    if abs(tr_l) > 1e-300:
        l_AR = l_AR / tr_l
    return r_AL, l_AR


def environments(AL_list, AR_list, W_list, D):
    """`(GL, GR, e_cell)` -- the per-bond channel-resolved environments for
    the whole cell, and the energy density per unit cell.

    `GL[n]` is the environment on the bond immediately to the LEFT of site
    `n` (so `GL[0]` sits at the cell's left edge), `GR[n]` the one
    immediately to its RIGHT. Both are `(Dw, D, D)`. They are built once at
    the cell edges and then simply pushed across, which is the whole point
    of the sequential formulation: one cell traversal gives every site's
    environment, rather than one solve per site."""
    Dw = W_list[0].shape[0]
    n_uc = len(W_list)
    r_AL, l_AR = _cell_fixed_points(AL_list, AR_list, D)

    GL0, src_l = _nilpotent_channels_left(AL_list, W_list, D, Dw)
    GRlast, src_r = _nilpotent_channels_right(AR_list, W_list, D, Dw)

    e_L = float(np.trace(r_AL.conj() @ src_l).real)
    e_R = float(np.trace(l_AR.conj() @ src_r).real)
    GL0 = GL0.copy()
    GL0[_F_IDX] = _solve_left_environment_cell(AL_list, r_AL, src_l, e_L, D)
    GRlast = GRlast.copy()
    GRlast[_S_IDX] = _solve_right_environment_cell(AR_list, l_AR, src_r, e_R, D)

    GL = [GL0]
    for n in range(n_uc - 1):
        GL.append(_push_cell_left(GL[n], [AL_list[n]], [W_list[n]]))
    GR = [None] * n_uc
    GR[n_uc - 1] = GRlast
    for n in range(n_uc - 1, 0, -1):
        GR[n - 1] = _push_cell_right(GR[n], [AR_list[n]], [W_list[n]])
    return GL, GR, 0.5 * (e_L + e_R)


# -- effective Hamiltonians --------------------------------------------------

def h_ac_action(X, GLn, GRn, Wn):
    """`H_AC[n]` applied to `X` (D,d,D):

        Y[L,o,R] = sum GL[a][l,L] X[l,i,r] W[a,i,o,b] GR[b][r,R]

    the standard one-site effective Hamiltonian, with both environments
    channel-resolved so every term (on-site, the two half-open bond terms,
    and the accumulated background on either side) is a single contraction
    rather than the four separate diagrams the reach-1 special case needs."""
    Dw = Wn.shape[0]
    Y = np.zeros_like(X)
    for a in range(Dw):
        if not GLn[a].any():
            continue
        # X with its left leg closed against GL[a]: einsum('lL,lir->Lir')
        XL = np.tensordot(GLn[a], X, axes=([0], [0]))
        for b in range(Dw):
            slab = Wn[a, :, :, b]
            if not slab.any() or not GRn[b].any():
                continue
            XLW = idmrg_exc._apply_op_ket(slab, XL)
            Y += np.tensordot(XLW, GRn[b], axes=([2], [0]))
    return Y


def h_c_action(C, GL_right_of_n, GR_n):
    """`H_C[n]` applied to the bond matrix `C` (D,D):

        Y[L,R] = sum_a GL[a][l,L] C[l,r] GR[a][r,R]

    the zero-site effective Hamiltonian. Both environments are taken ON
    that bond -- the left one is the environment to the left of site n+1,
    the right one is the environment to the right of site n -- and the same
    channel index closes on both sides, which is what makes this a plain
    sum of outer products."""
    Dw = GL_right_of_n.shape[0]
    Y = np.zeros_like(C)
    for a in range(Dw):
        L, R = GL_right_of_n[a], GR_n[a]
        if not L.any() or not R.any():
            continue
        Y += L.T @ C @ R
    return Y


# -- gauge update ------------------------------------------------------------

def update_AL_AR(AC, C_left, C_right):
    """New `(AL, AR)` for one site from its `AC` and the bond matrices on
    either side, via the orthogonal-Procrustes update: `AL` minimizes
    `||AC - AL @ C_right||` over isometric `AL`, `AR` minimizes
    `||AC - C_left @ AR||`.

    The multi-site case is where `vumps.py`'s single-site version needs
    generalizing rather than merely repeating: there, one `C` sits on both
    sides of the only site, so the same polar factor serves twice. Here the
    two bonds are different matrices, and using one for both silently
    gauges the state wrong."""
    D, d, _ = AC.shape
    U_l, _ = _polar(AC.reshape(D * d, D), side='right')
    U_cr, _ = _polar(C_right, side='right')
    AL_new = (U_l @ U_cr.conj().T).reshape(D, d, D)

    U_r, _ = _polar(AC.reshape(D, d * D), side='left')
    U_cl, _ = _polar(C_left, side='right')
    AR_new = (U_cl.conj().T @ U_r).reshape(D, d, D)
    return AL_new, AR_new


def gauge_mismatch(AC_list, C_list, AL_list, AR_list):
    """`max_n (||AC[n] - AL[n]C[n]|| + ||AC[n] - C[n-1]AR[n]||)/||AC[n]||`
    -- zero exactly at a fixed point. The max over the cell, not the mean:
    one badly gauged site in the cell is not converged, however well the
    others behave."""
    n_uc = len(AC_list)
    worst = 0.0
    for n in range(n_uc):
        D, d, _ = AC_list[n].shape
        left = AL_list[n].reshape(D * d, D) @ C_list[n]
        right = C_list[(n - 1) % n_uc] @ AR_list[n].reshape(D, d * D)
        nrm = np.linalg.norm(AC_list[n])
        if nrm <= 0:
            continue
        m = (np.linalg.norm(AC_list[n].reshape(D * d, D) - left)
             + np.linalg.norm(AC_list[n].reshape(D, d * D) - right)) / nrm
        worst = max(worst, float(m))
    return worst


# -- initial state -----------------------------------------------------------

def random_initial_state(D, dims, rng=None):
    """`(AL, AR, C)` lists for a cell whose sites have physical dimensions
    `dims`. Each `AL[n]` is made left-canonical and each `AR[n]`
    right-canonical individually, exactly as `vumps.py`'s own
    `_random_initial_state` does for its single supersite -- the two are
    deliberately NOT mutually consistent, which is the whole premise of the
    VUMPS iteration (it is the gauge mismatch between them that drives the
    updates, and starting from a consistent pair would be starting from a
    fixed point of nothing)."""
    if rng is None:
        rng = np.random.default_rng()
    AL, AR, C = [], [], []
    for d in dims:
        M = rng.standard_normal((D * d, D)) + 1j * rng.standard_normal((D * d, D))
        Q, _ = np.linalg.qr(M)
        AL.append(Q[:, :D].reshape(D, d, D))
        M2 = rng.standard_normal((D, d * D)) + 1j * rng.standard_normal((D, d * D))
        Q2, _ = np.linalg.qr(M2.conj().T)
        AR.append(Q2[:, :D].conj().T.reshape(D, d, D))
        C.append(np.eye(D, dtype=complex) / np.sqrt(D))
    return AL, AR, C


def grow_initial_state(D, dims, AL_old, AR_old, C_old, rng=None):
    """The same, warm-started from a converged smaller-`D` solution: embed
    it in the top-left block and fill the rest with small noise, then
    re-canonicalize. Mirrors `vumps.py`'s `_grow_initial_state`, and exists
    for the same measured reason -- a pure random start at large `D` lands
    in a bad basin often enough to matter, so the driver ramps `D` and
    warm-starts each step from the last."""
    if rng is None:
        rng = np.random.default_rng()
    D_old = AL_old[0].shape[0]
    AL, AR, C = [], [], []
    for n, d in enumerate(dims):
        def _embed(T_old):
            T = 0.05 * (rng.standard_normal((D, d, D))
                        + 1j * rng.standard_normal((D, d, D)))
            k = min(D, D_old)
            T[:k, :, :k] = T_old[:k, :, :k]
            return T
        A = _embed(AL_old[n])
        Q, _ = np.linalg.qr(A.reshape(D * d, D))
        AL.append(Q[:, :D].reshape(D, d, D))
        B = _embed(AR_old[n])
        Q2, _ = np.linalg.qr(B.reshape(D, d * D).conj().T)
        AR.append(Q2[:, :D].conj().T.reshape(D, d, D))
        Cn = np.zeros((D, D), dtype=complex)
        k = min(D, D_old)
        Cn[:k, :k] = C_old[n][:k, :k]
        Cn += 1e-3 * (rng.standard_normal((D, D)) + 1j * rng.standard_normal((D, D)))
        C.append(Cn / max(np.linalg.norm(Cn), 1e-300))
    return AL, AR, C


# -- the iteration -----------------------------------------------------------

def single_run(W_list, dims, D, tol, maxiter, niter_lanczos, init=None,
                rng=None):
    """One sequential multi-site VUMPS attempt.

    Each iteration builds the environments once for the whole cell and then
    sweeps it, solving `H_AC[n]` at every site and `H_C[n]` at every bond,
    and updating `(AL[n], AR[n])` from the results. That is the "sequential"
    structure the linear-scaling literature describes (arXiv:2003.01142) and
    that TeNPy implements as a `Sweep` over the unit cell -- as opposed to
    one grouped eigenproblem of dimension `D * prod(d_p) * D`.

    Returns a dict with the converged tensors, environments and diagnostics
    (a plain dict rather than a class: `vumps.py` owns the public
    `VUMPSResult`, and this is the piece it wraps)."""
    if rng is None:
        rng = np.random.default_rng()
    n_uc = len(W_list)
    if init is None:
        AL, AR, C = random_initial_state(D, dims, rng)
    else:
        AL, AR, C = [list(x) for x in init]

    AC = [AL[n].reshape(-1, D) @ C[n] for n in range(n_uc)]
    AC = [AC[n].reshape(D, dims[n], D) for n in range(n_uc)]

    converged, mismatch, it = False, float('inf'), 0
    GL = GR = None
    e_cell = 0.0
    for it in range(maxiter):
        GL, GR, e_cell = environments(AL, AR, W_list, D)

        AC_new, C_new = [None] * n_uc, [None] * n_uc
        for n in range(n_uc):
            def _hac(x, n=n):
                return h_ac_action(x.reshape(D, dims[n], D), GL[n], GR[n],
                                    W_list[n]).reshape(-1)
            _e, v = _lanczos_ground_state(_hac, AC[n].reshape(-1),
                                           niter=niter_lanczos)
            AC_new[n] = v.reshape(D, dims[n], D)
            # H_C[n] lives on the bond to the right of site n: its left
            # environment is the one to the left of site n+1 (which is
            # GL[n+1], or the next cell's GL[0] when n is the last site),
            # its right environment is GR[n].
            GL_bond = GL[(n + 1) % n_uc] if n + 1 < n_uc else \
                _push_cell_left(GL[n], [AL[n]], [W_list[n]])

            def _hc(x, GL_bond=GL_bond, n=n):
                return h_c_action(x.reshape(D, D), GL_bond, GR[n]).reshape(-1)
            _e2, vc = _lanczos_ground_state(_hc, C[n].reshape(-1),
                                             niter=niter_lanczos)
            C_new[n] = vc.reshape(D, D)

        mismatch = gauge_mismatch(AC_new, C_new, AL, AR)
        AC, C = AC_new, C_new
        for n in range(n_uc):
            AL[n], AR[n] = update_AL_AR(AC[n], C[(n - 1) % n_uc], C[n])
        if mismatch < tol:
            converged = True
            break

    # Refresh against the FINAL AL/AR: the environments and energy above
    # were built from this iteration's INPUT tensors, one update behind
    # what is being returned. Same reason vumps.py's own single run does a
    # final `_environments` call after its loop.
    GL, GR, e_cell = environments(AL, AR, W_list, D)
    return dict(AL=AL, AR=AR, C=C, AC=AC, GL=GL, GR=GR, e_cell=e_cell,
                converged=converged, niter=min(it + 1, maxiter),
                mismatch=mismatch)


def _d_ramp(D):
    """1, 2, 4, ..., D -- the bond dimensions solved at on the way to `D`.
    Same doubling ladder as `vumps.py`'s own `_d_ramp`, and for the same
    reason: a pure random start at large `D` lands in a bad basin often
    enough to matter, while ramping one integer at a time would make the
    driver O(D^4) with the ramp, not the target solve, dominating."""
    ramp, d = [], 1
    while d < D:
        ramp.append(d)
        d *= 2
    ramp.append(D)
    return ramp


def ground_state(W_list, dims, D, tol=1e-10, maxiter=800, niter_lanczos=40,
                  nrestarts=4, verbose=False, rng=None):
    """Sequential multi-site VUMPS with the same restart/D-ramp discipline
    `vumps.py`'s own driver uses -- warm-start each rung of the ramp from
    the previous one, try `nrestarts` attempts per rung, keep the best
    (converged beats non-converged, then lowest energy), and spend a
    bounded extra budget if a larger `D` came out ABOVE an already-known
    smaller-`D` energy, which the variational principle forbids and which
    is a real, reachable outcome of a restart search rather than a
    hypothetical.

    Returns the same dict `single_run` does."""
    if rng is None:
        rng = np.random.default_rng()

    def attempt(D_cur, init):
        try:
            return single_run(W_list, dims, D_cur, tol, maxiter,
                               niter_lanczos, init=init, rng=rng)
        except Exception as exc:                      # noqa: BLE001
            # A degenerate transfer-matrix spectrum or a singular
            # regularized solve: recoverable by retrying from a different
            # random start, exactly as vumps.py's own driver treats it.
            if verbose:
                print("vumps_ms D={} attempt failed: {!r}".format(D_cur, exc))
            return None

    def better(a, b):
        if a is None:
            return b
        if b is None:
            return a
        if a["converged"] != b["converged"]:
            return a if a["converged"] else b
        return a if a["e_cell"] < b["e_cell"] else b

    best, prev, best_e = None, None, None
    for D_cur in _d_ramp(D):
        n_here = nrestarts if D_cur == D else min(nrestarts, 3)
        local = None
        for i in range(n_here):
            init = None
            if i == 0 and prev is not None:
                init = grow_initial_state(D_cur, dims, prev["AL"], prev["AR"],
                                           prev["C"], rng)
            local = better(attempt(D_cur, init), local)
        if local is None:
            raise RuntimeError(
                "vumps_ms.ground_state: every attempt at D={} failed -- try "
                "increasing nrestarts".format(D_cur))
        if best_e is not None and local["e_cell"] > best_e + 1e-6:
            for _ in range(2 * nrestarts):
                if local["e_cell"] <= best_e + 1e-6:
                    break
                init = (grow_initial_state(D_cur, dims, prev["AL"], prev["AR"],
                                            prev["C"], rng)
                        if prev is not None else None)
                local = better(attempt(D_cur, init), local)
        if verbose:
            print("vumps_ms D={}: e_cell={} converged={}".format(
                D_cur, local["e_cell"], local["converged"]))
        prev = local
        best_e = local["e_cell"] if best_e is None else min(best_e, local["e_cell"])
        best = local
    return best


# -- observables -------------------------------------------------------------

def onsite_expectation(AC_list, sites_uc, opname, p):
    """`<opname>` at site `p` of the cell.

    `AC[p]` is already the exactly normalized single-site reduced state by
    construction of the mixed canonical gauge (everything to its left
    contracts to the identity through `AL`, everything to its right through
    `AR`) -- Vanderstraeten, Haegeman & Verstraete, arXiv:1810.07006
    Eq.(34) -- so this is a plain contraction, with no eigenproblem and no
    grouping. Same formula `vumps.onsite_expectation` uses; the only
    difference is that the relevant `AC` is picked out of the cell rather
    than embedded into a supersite operator."""
    M = sites_uc.site_type(p + 1).matrix(opname)
    AC = AC_list[p]
    AC_op = idmrg_exc._apply_op_ket(M, AC)
    val = np.tensordot(np.conj(AC), AC_op, axes=([0, 1, 2], [0, 1, 2]))
    nrm = np.tensordot(np.conj(AC), AC, axes=([0, 1, 2], [0, 1, 2])).real
    return complex(val / nrm)


def two_point_correlator(AC_list, AR_list, sites_uc, n_uc, opname_i, p_i,
                          opname_j, r):
    """`<opname_i(p_i) opname_j(p_i + r)>`, `r` in physical sites.

    Puts the mixed-gauge centre at `p_i` (so `AC[p_i]` carries the exact
    normalization and everything to its left is already the identity), then
    walks right through `AR` -- which is exactly right-orthonormal, so the
    tail beyond the second operator closes to the identity with no
    eigenproblem either. Positions advance through the cell modulo `n_uc`,
    and the site TYPE at each position is set by `position % n_uc`.

    The endpoint matrices, the ordering sign and whether a Jordan-Wigner
    string is open all come from `idmrg._term_site_matrices` -- the same
    helper that builds the Hamiltonian's own two-site terms -- so a
    correlator and the Hamiltonian cannot disagree about the fermionic
    convention. Parity-even operators come back stringless, reproducing the
    plain bosonic formula exactly."""
    if r < 0:
        raise ValueError("two_point_correlator: r must be >= 0, got {}".format(r))
    term = [1.0, [opname_i, p_i], [opname_j, p_i + r]]
    _rel, coef, mats, ferm = idmrg._term_site_matrices(term, sites_uc, n_uc)
    strung = ferm[0] if r > 0 else False

    if r == 0:
        # Same-site: compose the two matrices (M_j @ M_i, idmrg's own r=0
        # convention) and measure as a one-site expectation.
        M = mats[1] @ mats[0] if len(mats) > 1 else mats[0]
        AC = AC_list[p_i % n_uc]
        AC_op = idmrg_exc._apply_op_ket(M, AC)
        val = np.tensordot(np.conj(AC), AC_op, axes=([0, 1, 2], [0, 1, 2]))
        nrm = np.tensordot(np.conj(AC), AC, axes=([0, 1, 2], [0, 1, 2])).real
        return complex(coef * val / nrm)

    AC = AC_list[p_i % n_uc]
    X = _push_left(np.eye(AC.shape[0], dtype=complex), AC, mats[0])
    # Divide out AC's own norm once, so the walk below closes at 1.
    nrm = np.tensordot(np.conj(AC), AC, axes=([0, 1, 2], [0, 1, 2])).real
    for step in range(1, r + 1):
        pos = p_i + step
        A = AR_list[pos % n_uc]
        if step == r:
            M = mats[1]
        elif strung:
            M = sites_uc.site_type(pos % n_uc + 1).matrix("F")
        else:
            M = None
        X = _push_left(X, A, M)
    return complex(coef * np.trace(X) / nrm)
