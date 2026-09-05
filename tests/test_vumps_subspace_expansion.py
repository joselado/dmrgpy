"""Coverage for pyitensor/vumps.py's subspace expansion -- the way the
D-ramp grows an already-converged solution at D_old into a starting point
at D_new (`vumps._subspace_expand`, modelled on ITensorInfiniteMPS.jl's
own `subspace_expansion.jl`; see that function's own section comment).

What is pinned here is the set of invariants that make the expansion safe
to put in front of the ramp at all, and every one of them is exact (not
approximate), so these assertions are tight rather than tolerance-tuned:

 1. the enlarged AL/AR are still exactly isometric (left/right canonical);
 2. the enlarged state is the SAME physical state -- AC is unchanged in
    the old block and exactly zero in the new one, so the expansion cannot
    move the energy;
 3. the energy density computed from the enlarged tensors therefore equals
    the pre-expansion one to machine precision;
 4. the new directions are genuinely orthogonal to the old ones (this is
    what item 1 rests on, checked directly rather than only through it);
 5. the ramp still reaches the right answer end to end.

Item 3 is the one that matters physically: a warm start that changed the
energy would break the variational-principle safety net in
`vumps_ground_state` (which compares each ramp step against the best
energy already achieved at a smaller D).
"""
import numpy as np
import pytest
from scipy.integrate import quad

from dmrgpy import infinitechain
from dmrgpy.pyitensor import idmrg
from dmrgpy.pyitensor import idmrg_excitations as idmrg_exc
from dmrgpy.pyitensor import vumps


def _heisenberg():
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
    ic.set_hamiltonian(ic.SxC[0] * ic.SxR[0] + ic.SyC[0] * ic.SyR[0]
                       + ic.SzC[0] * ic.SzR[0])
    return ic


def _tfim(g):
    """H = -sum(sigma^x sigma^x) - g*sum(sigma^z), as in tests/test_vumps.py."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
    ic.set_hamiltonian(-4.0 * ic.SxC[0] * ic.SxR[0] - 2.0 * g * ic.SzC[0])
    return ic


def _tfim_exact_energy_density(g):
    val, _ = quad(lambda k: np.sqrt(1 + g ** 2 - 2 * g * np.cos(k)), 0, np.pi)
    return -val / np.pi


def _solve_at(ic, D, seed=0):
    """One raw `_vumps_single_run` at bond dimension D, plus the grouped
    automaton pieces `_subspace_expand` itself needs."""
    sites_uc, W_bulk = idmrg._build_automaton(
        ic._h_intra.op, ic._h_inter.op, ic.site_types, ic.n_uc)
    W = vumps._group_automaton(W_bulk, ic.n_uc)
    idmrg_exc._check_reach_one(W)
    d_g = int(np.prod([sites_uc.dim(p + 1) for p in range(ic.n_uc)]))
    pending = idmrg_exc._pending_channels(W)
    h1 = idmrg_exc._onsite_matrix(W)
    np.random.seed(seed)
    res = vumps._vumps_single_run(sites_uc, ic.n_uc, D, d_g, W, pending, h1,
                                  1e-10, 800, 30, False)
    return res, W, pending, h1


@pytest.mark.parametrize("D_old,D_new", [(1, 2), (2, 4)])
def test_expansion_preserves_isometry_and_state(D_old, D_new):
    ic = _tfim(1.5)
    res, W, pending, h1 = _solve_at(ic, D_old)
    expanded = vumps._subspace_expand(res, pending, h1, D_new)
    assert expanded is not None, "expansion found no direction to grow into"
    AL, AR, C = expanded
    Dn, d_g, _ = AL.shape
    assert Dn == D_new

    # (1) still exactly canonical
    left_gram = np.einsum('lpm,lpn->mn', AL.conj(), AL)
    right_gram = np.einsum('lpr,mpr->lm', AR, AR.conj())
    assert np.abs(left_gram - np.eye(Dn)).max() == pytest.approx(0.0, abs=1e-12)
    assert np.abs(right_gram - np.eye(Dn)).max() == pytest.approx(0.0, abs=1e-12)

    # (4) the added columns/rows are the ones the old isometries did not reach
    assert np.abs(np.einsum('lpm,lpa->ma', res.AL.conj(),
                            AL[:D_old, :, D_old:])).max() == pytest.approx(0.0, abs=1e-12)
    assert np.abs(np.einsum('apr,lpr->al', AR[D_old:, :, :D_old],
                            res.AR.conj())).max() == pytest.approx(0.0, abs=1e-12)

    # (2) same physical state: AC unchanged where it existed, zero elsewhere
    AC_old = np.einsum('lpm,mn->lpn', res.AL, res.C)
    AC_new = np.einsum('lpm,mn->lpn', AL, C)
    assert np.abs(AC_new[:D_old, :, :D_old] - AC_old).max() == pytest.approx(0.0, abs=1e-14)
    assert np.abs(AC_new[:, :, D_old:]).max() == pytest.approx(0.0, abs=1e-14)
    assert np.abs(AC_new[D_old:, :, :]).max() == pytest.approx(0.0, abs=1e-14)

    # (3) and therefore the same energy density
    _GL, _GR, e_cell, _bond_envs = vumps._environments(AL, AR, W, pending)
    assert e_cell / ic.n_uc == pytest.approx(res.e0, abs=1e-10)


def test_expansion_declines_when_there_is_nothing_to_grow_into():
    """D_new <= D is not an expansion, and must be reported as such (the
    driver falls back to the noise start on a None) rather than returning
    a same-size or shrunken state."""
    ic = _tfim(1.5)
    res, _W, pending, h1 = _solve_at(ic, 2)
    assert vumps._subspace_expand(res, pending, h1, 2) is None
    assert vumps._subspace_expand(res, pending, h1, 1) is None


def test_expansion_is_capped_by_the_null_space_size():
    """d_g=2 leaves only D*(d_g-1)=D new directions, so a ramp step asking
    for more than a doubling gets back a SMALLER state than it asked for --
    the caller (`vumps_ground_state.warm_start`) is what pads the rest, so
    this must not silently over-report its own dimension."""
    ic = _tfim(1.5)
    res, _W, pending, h1 = _solve_at(ic, 2)
    AL, _AR, _C = vumps._subspace_expand(res, pending, h1, 16)
    assert AL.shape[0] == 4


def test_ramp_with_expansion_reaches_the_exact_tfim_energy():
    """End to end through the public driver, whose ramp now warm-starts
    every step by expansion: a gapped TFIM at D=4 should sit close to the
    exact free-fermion energy density (Pfeuty 1970)."""
    ic = _tfim(1.5)
    np.random.seed(0)
    res = vumps.vumps_ground_state(ic.site_types, ic._h_intra.op,
                                   ic._h_inter.op, ic.n_uc, 4,
                                   tol=1e-10, maxiter=800, nrestarts=3)
    assert res.e0 == pytest.approx(_tfim_exact_energy_density(1.5), abs=1e-5)


@pytest.mark.parametrize("model", ["tfim", "heisenberg"])
def test_two_site_effective_hamiltonian_is_hermitian(model):
    """`vumps._h_two_site_action` assembled densely over a basis must be
    Hermitian -- H is, and every piece of the two-site widening (both
    onsite terms, the bond term inside the pair, and the two bond terms
    that straddle the pair's edges into the precomputed neighbour
    environments) preserves that.

    Heisenberg is the case that has to be here. Every other assertion in
    this file -- isometry, state preservation, energy preservation -- holds
    whatever operator the direction heuristic is built from, so none of
    them can see a transposed factor in the intra-pair term; and TFIM
    cannot either, because its own bond operators are the symmetric Sx.
    Heisenberg's Sy is antisymmetric, so it fails immediately (this is not
    hypothetical -- it is exactly the bug this test was written for)."""
    ic = _heisenberg() if model == "heisenberg" else _tfim(1.5)
    res, _W, pending, h1 = _solve_at(ic, 3)
    bond_envs = vumps._precompute_bond_environments(res.AL, res.AR, pending)
    D, d_g, _ = res.AL.shape
    n = D * d_g * d_g * D
    M = np.zeros((n, n), dtype=complex)
    for k in range(n):
        v = np.zeros(n, dtype=complex)
        v[k] = 1.0
        M[:, k] = vumps._h_two_site_action(
            v.reshape(D, d_g, d_g, D), res.GL, res.GR, bond_envs, h1).reshape(-1)
    scale = max(np.abs(M).max(), 1.0)
    assert np.abs(M - M.conj().T).max() / scale == pytest.approx(0.0, abs=1e-12)


def test_undersized_expansion_is_paddable_to_the_requested_dimension():
    """The fallback path in `vumps_ground_state`'s own `warm_start`: when
    the null space cannot supply the whole ramp step, the (smaller)
    expanded tensors are handed to `_grow_initial_state`, which noise-pads
    the rest of the way. It has to accept them -- they are exactly what a
    D_old<D seed looks like to it."""
    ic = _tfim(1.5)
    res, _W, pending, h1 = _solve_at(ic, 2)
    AL_e, AR_e, _C_e = vumps._subspace_expand(res, pending, h1, 6)
    assert AL_e.shape[0] == 4 < 6
    AL, AR, C = vumps._grow_initial_state(6, AL_e.shape[1], AL_e, AR_e)
    assert AL.shape == (6, 2, 6) and AR.shape == (6, 2, 6) and C.shape == (6, 6)
