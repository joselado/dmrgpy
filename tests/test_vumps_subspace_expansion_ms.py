"""Coverage for pyitensor/vumps_ms.py's PER-BOND subspace expansion -- the
multi-site counterpart of `vumps._subspace_expand`, and the shape
ITensorInfiniteMPS.jl's own `subspace_expansion(psi, H)` has (a loop over
every bond of the unit cell, feeding an outer expand-then-solve ramp).

`vumps.py`'s version only ever runs at n_uc <= 2, where the cell is grouped
into one supersite and so has exactly one bond; ungrouped, a cell has n_uc
of them and the two null spaces at a bond come from DIFFERENT tensors. What
is checked here is therefore both that the generalization is right and that
it degenerates to the already-validated one:

 1. the new two-site effective Hamiltonian (`vumps_ms.h_two_site_action`,
    built from the channel-resolved environments) reproduces the grouped
    module's own `_h_two_site_action` (built from its reach-1 {GL, GR,
    bond_envs} triple) to machine precision at n_uc=1 -- on the physical
    theta AND on a random tensor, so it is the operator that agrees, not
    just its value on one vector;
 2. the expansion picks the SAME directions the grouped one does at n_uc=1
    (compared as a subspace, since the SVD fixes each direction only up to
    a phase);
 3. every invariant `tests/test_vumps_subspace_expansion.py` pins for the
    grouped version holds per site here: the enlarged AL/AR stay exactly
    isometric, AC is unchanged, and the energy density does not move;
 4. the ramp still reaches the right answer end to end at n_uc >= 3.

Item 3 is the one that matters physically: a warm start that moved the
energy would break `ground_state`'s variational-principle safety net, which
compares each ramp step against the best energy already reached.
"""
import numpy as np
import pytest
from scipy.integrate import quad

from dmrgpy.pyitensor import idmrg
from dmrgpy.pyitensor import idmrg_excitations as idmrg_exc
from dmrgpy.pyitensor import vumps
from dmrgpy.pyitensor import vumps_ms

SPIN_HALF = 2          # sites/siteset.py's own type code


def _heisenberg_terms(n_uc, J=1.0):
    intra, inter = [], []
    for i in range(n_uc):
        for op in ("Sx", "Sy", "Sz"):
            term = [J, [op, i], [op, i + 1]]
            (intra if i + 1 < n_uc else inter).append(term)
    return intra, inter


def _tfim_terms(n_uc, g=1.5):
    intra, inter = [], []
    for i in range(n_uc):
        term = [-4.0, ["Sx", i], ["Sx", i + 1]]
        (intra if i + 1 < n_uc else inter).append(term)
        intra.append([-2.0 * g, ["Sz", i]])
    return intra, inter


def _tfim_exact_energy_density(g):
    val, _ = quad(lambda k: np.sqrt(1 + g ** 2 - 2 * g * np.cos(k)), 0, np.pi)
    return -val / np.pi


def _automaton(terms, n_uc):
    intra, inter = terms(n_uc)
    sites_uc, W_bulk = idmrg._build_automaton(
        intra, inter, [SPIN_HALF] * n_uc, n_uc)
    W_list = [W_bulk[p].array for p in range(n_uc)]
    dims = [sites_uc.dim(p + 1) for p in range(n_uc)]
    return sites_uc, W_bulk, W_list, dims


def _solve_at(terms, n_uc, D, seed=0, maxiter=300):
    sites_uc, W_bulk, W_list, dims = _automaton(terms, n_uc)
    out = vumps_ms.single_run(W_list, dims, D, 1e-10, maxiter, 40,
                               rng=np.random.default_rng(seed))
    return out, W_list, dims, sites_uc, W_bulk


@pytest.mark.parametrize("terms", [_heisenberg_terms, _tfim_terms])
def test_two_site_action_matches_the_grouped_one(terms):
    """The channel-resolved two-site H must BE the grouped module's own
    four-diagram one at n_uc=1, not merely agree on the ground state: the
    random-tensor case is what makes this an operator identity."""
    out, W_list, dims, sites_uc, W_bulk = _solve_at(terms, 1, 3, seed=7)
    Wg = vumps._group_automaton(W_bulk, 1)
    pending = idmrg_exc._pending_channels(Wg)
    h1 = idmrg_exc._onsite_matrix(Wg)
    GLg, GRg, _e, bond_envs = vumps._environments(
        out["AL"][0], out["AR"][0], Wg, pending)

    D, d = out["AL"][0].shape[0], dims[0]
    theta = np.einsum('lpx,xy,yqr->lpqr', out["AL"][0], out["C"][0], out["AR"][0])
    rng = np.random.default_rng(3)
    random_theta = (rng.standard_normal((D, d, d, D))
                    + 1j * rng.standard_normal((D, d, d, D)))
    for th in (theta, random_theta):
        got = vumps_ms.h_two_site_action(th, out["GL"][0], out["GR"][0],
                                          W_list[0], W_list[0])
        want = vumps._h_two_site_action(th, GLg, GRg, bond_envs, h1)
        assert np.max(np.abs(got - want)) < 1e-12 * max(np.max(np.abs(want)), 1.0)


def test_expansion_picks_the_same_directions_as_the_grouped_one():
    """At n_uc=1 both modules expand the same single bond. The SVD fixes
    each new direction only up to a phase, so what must agree is the
    SUBSPACE the added directions span, not the tensors themselves.

    Note what this does and does not reach: at `D=4 -> 8` on a d=2 chain,
    `keep` equals the whole null-space dimension `D*(d-1)`, so both
    projectors are the identity on that space and it is the two NULL
    SPACES agreeing that is checked, not the SVD's selection among them.
    That is deliberate -- the selection is pinned instead by
    `test_the_cell_boundary_bond_uses_the_right_environments` below, which
    compares the singular values themselves."""
    D, D_new = 4, 8
    out, W_list, dims, sites_uc, W_bulk = _solve_at(_heisenberg_terms, 1, D, seed=5)
    Wg = vumps._group_automaton(W_bulk, 1)
    pending = idmrg_exc._pending_channels(Wg)
    h1 = idmrg_exc._onsite_matrix(Wg)
    GLg, GRg, e_g, bond_envs = vumps._environments(
        out["AL"][0], out["AR"][0], Wg, pending)
    grouped = vumps.VUMPSResult(sites_uc, 1, D, dims[0], out["AL"][0],
                                 out["AR"][0], out["C"][0], out["AC"][0],
                                 GLg, GRg, Wg, e_g, True, 1, 0.0)

    AL_g, AR_g, _C_g = vumps._subspace_expand(grouped, pending, h1, D_new)
    AL_m, AR_m, _C_m = vumps_ms.subspace_expand(out, W_list, dims, D_new)
    assert AL_g.shape == AL_m[0].shape

    k = AL_g.shape[0] - D
    add_g = AL_g[:D, :, D:].reshape(-1, k)
    add_m = AL_m[0][:D, :, D:].reshape(-1, k)
    assert np.max(np.abs(add_g @ add_g.conj().T
                         - add_m @ add_m.conj().T)) < 1e-10
    add_g = AR_g[D:, :, :D].reshape(k, -1)
    add_m = AR_m[0][D:, :, :D].reshape(k, -1)
    assert np.max(np.abs(add_g.conj().T @ add_g
                         - add_m.conj().T @ add_m)) < 1e-10


def test_the_cell_boundary_bond_uses_the_right_environments():
    """The cell's LAST bond is the one piece of the multi-site
    generalization the n_uc=1 oracle above cannot reach: it straddles the
    cell boundary, so its right environment/automaton/tensor are site 0's
    of the NEXT cell, taken to be this cell's own `GR[0]`/`W_list[0]`/
    `AR[0]` by periodicity. Nothing else here would catch that being
    wrong -- the isometry/state/energy invariants hold whatever `H2`
    computes, since it only picks directions.

    What pins it: the singular values of `M = NL^dagger (H2 theta)
    NR^dagger` are gauge invariant (a bond gauge `C -> U C V^dagger` sends
    `M -> U M V^dagger`), and a constant shift of `H` -- the per-cell
    energy baseline the environments carry, which differs between a 1- and
    a 2-site cell -- contributes `c*theta`, whose overlap with the two null
    spaces is exactly zero. So on one uniform chain, the boundary bond of a
    2-site cell, its interior bond, and the single bond of a 1-site cell
    must all give the same singular values."""
    D = 4
    per_cell = {}
    for n_uc in (1, 2):
        _sites, _W_bulk, W_list, dims = _automaton(_tfim_terms, n_uc)
        out = vumps_ms.ground_state(W_list, dims, D, tol=1e-12, maxiter=600,
                                     niter_lanczos=40, nrestarts=3,
                                     rng=np.random.default_rng(0))
        assert out["converged"]
        svs = []
        for n in range(n_uc):
            m = (n + 1) % n_uc
            NL = vumps._null_space_left(out["AL"][n])
            NR = vumps._null_space_right(out["AR"][m])
            theta = np.einsum('lpx,xy,yqr->lpqr', out["AL"][n], out["C"][n],
                              out["AR"][m])
            H2theta = vumps_ms.h_two_site_action(theta, out["GL"][n],
                                                  out["GR"][m], W_list[n],
                                                  W_list[m])
            M = np.einsum('lpa,lpqr,bqr->ab', NL.conj(), H2theta, NR.conj())
            svs.append(np.linalg.svd(M, compute_uv=False))
        per_cell[n_uc] = (out["e_cell"] / n_uc, svs)

    (e1, svs1), (e2, svs2) = per_cell[1], per_cell[2]
    assert e1 == pytest.approx(e2, abs=1e-10)      # the same state, first
    k = min(len(svs1[0]), len(svs2[0]))
    for bond, sv in enumerate(svs2):               # bond 1 is the boundary one
        assert np.max(np.abs(sv[:k] - svs1[0][:k])) < 1e-10, \
            "n_uc=2 bond {} disagrees with the 1-site cell".format(bond)


@pytest.mark.parametrize("n_uc", [1, 2, 3, 4])
def test_expansion_preserves_isometry_and_state(n_uc):
    """Exact invariants, so these are machine-precision assertions: the
    enlarged tensors are still canonical, the state is untouched (AC is
    the old one embedded), and the energy density therefore cannot move."""
    D, D_new = 3, 6
    out, W_list, dims, _sites, _W_bulk = _solve_at(_heisenberg_terms, n_uc, D,
                                                    seed=2, maxiter=200)
    expanded = vumps_ms.subspace_expand(out, W_list, dims, D_new)
    assert expanded is not None
    AL, AR, C = expanded
    Dn = AL[0].shape[0]
    assert D < Dn <= D_new

    for n in range(n_uc):
        M = AL[n].reshape(-1, Dn)
        assert np.max(np.abs(M.conj().T @ M - np.eye(Dn))) < 1e-10
        R = AR[n].reshape(Dn, -1)
        assert np.max(np.abs(R @ R.conj().T - np.eye(Dn))) < 1e-10
        # AC = AL @ C: the old block is untouched and the new one is zero,
        # i.e. the added directions carry no weight at all.
        AC = (AL[n].reshape(-1, Dn) @ C[n]).reshape(Dn, dims[n], Dn)
        AC_old = (out["AL"][n].reshape(-1, D) @ out["C"][n]).reshape(D, dims[n], D)
        assert np.max(np.abs(AC[:D, :, :D] - AC_old)) < 1e-12
        assert np.max(np.abs(AC[D:, :, :])) < 1e-12
        assert np.max(np.abs(AC[:, :, D:])) < 1e-12

    _GL, _GR, e_new = vumps_ms.environments(AL, AR, W_list, Dn)
    assert e_new == pytest.approx(out["e_cell"], abs=1e-10)


def test_expansion_declines_when_there_is_nothing_to_grow_into():
    """A request that cannot grow the bond returns None rather than a
    same-sized state, which is the contract `_warm_start` relies on to
    fall back to noise-padding."""
    out, W_list, dims, _s, _w = _solve_at(_heisenberg_terms, 2, 4, seed=1,
                                           maxiter=100)
    assert vumps_ms.subspace_expand(out, W_list, dims, 4) is None
    assert vumps_ms.subspace_expand(out, W_list, dims, 2) is None


def test_warm_start_always_returns_the_requested_dimension():
    """`_warm_start` must hand `single_run` exactly D_cur whether the
    expansion filled it, partly filled it, or declined -- the module
    carries one D, not one per bond."""
    n_uc = 3
    out, W_list, dims, _s, _w = _solve_at(_heisenberg_terms, n_uc, 2, seed=4,
                                           maxiter=100)
    rng = np.random.default_rng(0)
    for D_cur in (3, 4, 8, 16):
        AL, AR, C = vumps_ms._warm_start(out, W_list, dims, D_cur, rng)
        assert len(AL) == len(AR) == len(C) == n_uc
        for n in range(n_uc):
            assert AL[n].shape == (D_cur, dims[n], D_cur)
            assert AR[n].shape == (D_cur, dims[n], D_cur)
            assert C[n].shape == (D_cur, D_cur)


def test_ramp_with_expansion_reaches_the_exact_tfim_energy():
    """End to end through the driver the expansion is actually wired into,
    on a model with a closed-form answer."""
    n_uc = 3
    _s, _w, W_list, dims = _automaton(_tfim_terms, n_uc)
    out = vumps_ms.ground_state(W_list, dims, 8, tol=1e-10, maxiter=400,
                                 niter_lanczos=40, nrestarts=2)
    assert out["e_cell"] / n_uc == pytest.approx(
        _tfim_exact_energy_density(1.5), abs=1e-6)
