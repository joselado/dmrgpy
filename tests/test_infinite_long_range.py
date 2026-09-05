"""Couplings reaching further than one unit cell.

`Infinite_Many_Body_Chain` used to reject any term touching both the
previous (L) and the next (R) cell -- "spans three cells, out of scope" --
so the longest coupling expressible on an n_uc-site cell reached n_uc sites.
A longer-range model therefore had to be rewritten on a cell at least as
long as its range, which `gs_method="vumps"` then folded into a d**n_uc
supersite: exponential in exactly the parameter that should be linear.

Nothing in the automaton ever needed that restriction --
`idmrg._build_periodic_mpo`/`_active_channels_at` have always carried one
pending channel per site of a term's reach -- and neither does the
sequential multi-site VUMPS solver (`pyitensor/vumps_ms.py`), whose
channel-resolved environments consume exactly those channels. So
`set_hamiltonian` now accepts any finite range, `get_operator(..., group=c)`
takes an integer cell offset to write one, and `vumps.vumps_ground_state`
routes such a Hamiltonian to the sequential solver at ANY n_uc.

The checks are of three kinds:

1. *Against an exact answer.* A strong field polarizes the chain, so the
   ground state is an exact product state at D=1 and every long-range
   Sz-Sz term contributes a known shift. A term silently dropped,
   double-counted or mis-tiled cannot survive this, at any reach.
2. *Cell-size invariance.* The same chain written on a 1-site cell
   (reach 2, the new sequential route) and on a 2-site cell (reach 1, the
   already-validated grouped route) must give the same energy density --
   two different code paths for one physical model.
3. *Canonicalization and gating.* That the term-position canonicalization
   is a translation (so no bond is double-counted or lost), and that the
   one algorithm which genuinely cannot do reach>1 -- the tangent-space
   excitation ansatz -- says so instead of answering.
"""
import numpy as np
import pytest

from dmrgpy import cppext
from dmrgpy import infinitechain
from dmrgpy.infinitechain import _canonicalize_hamiltonian, _window_hamiltonian

# The dispatch exists on both backends (vumps.vumps_ground_state's own
# reach test on the "python" side, Chain::vumps_ground_state's
# `use_multisite`/`vumps_is_reach_one` on the C++ one), so the physics
# checks that can run on both do.
BACKENDS = ["python"] + ([3] if cppext.available(3) else [])

FIELD, J2, J3 = 4.0, 0.7, 0.3
# Fully polarized (the field beats both couplings), so <Sz>=1/2 on every
# site: -FIELD*<Sz> + J2*<Sz><Sz> + J3*<Sz><Sz> per site.
POLARIZED_EXACT = -FIELD / 2.0 + J2 / 4.0 + J3 / 4.0


def _polarized_chain(n_uc, reaches, backend="python"):
    """-FIELD*sum Sz + sum_r J_r sum Sz_i Sz_{i+r}, on an n_uc-site cell."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"] * n_uc,
                                            itensor_version=backend)
    ic.maxm, ic.maxiter, ic.etol = 2, 300, 1e-12
    ic.vumps_nrestarts = 2
    h = 0
    for i in range(n_uc):
        h = h - FIELD * ic.SzC[i]
        for r, J in reaches:
            k = i + r
            other = (ic.SzC[k] if k < n_uc
                     else ic.get_operator("Sz", k % n_uc, group=k // n_uc))
            h = h + J * ic.SzC[i] * other
    ic.set_hamiltonian(h)
    return ic


def _ising_j2_chain(n_uc, g=2.5, J1=1.0, J2c=0.5, D=8, backend="python"):
    """H = -g sum sigma^z - J1 sum sigma^x_i sigma^x_{i+1}
           - J2 sum sigma^x_i sigma^x_{i+2}, deep in the paramagnetic phase.

    Gapped and uniform on purpose. A near-critical model (J1-J2 Heisenberg
    at J2/J1 ~ 0.4 was tried first) makes VUMPS's own non-convex restart
    search the thing under test rather than the reach machinery: it landed
    in different basins run to run and the comparison flaked, at ~100s a
    go. Here both cells converge to the same answer every time, in ~2s."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"] * n_uc,
                                            itensor_version=backend)
    ic.maxm, ic.maxiter, ic.etol = D, 300, 1e-10
    ic.vumps_nrestarts = 2
    h = 0
    for i in range(n_uc):
        h = h - 2.0 * g * ic.SzC[i]
        for r, J in ((1, J1), (2, J2c)):
            k = i + r
            other = (ic.SxC[k] if k < n_uc
                     else ic.get_operator("Sx", k % n_uc, group=k // n_uc))
            h = h - 4.0 * J * ic.SxC[i] * other
    ic.set_hamiltonian(h)
    return ic


# -- 1. exact answers --------------------------------------------------------

@pytest.mark.parametrize("n_uc", [1, 2, 3])
@pytest.mark.parametrize("gs_method", ["vumps", "idmrg"])
def test_polarized_chain_with_long_range_terms_is_exact(n_uc, gs_method):
    """Both ground-state algorithms handle reach>1, and this pins the
    coefficient of every long-range term rather than just its presence.
    (`idmrg` has its own, unrelated n_uc<=2 restriction.)"""
    if gs_method == "idmrg" and n_uc > 2:
        pytest.skip("idmrg_ground_state supports n_uc<=2 only")
    ic = _polarized_chain(n_uc, [(2, J2), (3, J3)])
    assert ic._reach_cells == max(1, -(-3 // n_uc))   # ceil(3/n_uc)
    ic.gs_method = gs_method
    assert ic.gs_energy() == pytest.approx(POLARIZED_EXACT, abs=1e-9)


@pytest.mark.parametrize("n_uc", [1, 2])
def test_static_observables_follow_the_long_range_dispatch(n_uc):
    """`vev`/`correlator` on a reach>1 chain -- which at n_uc<=2 means a
    `.multisite` VUMPSResult reaching `vumps_ms.onsite_expectation`/
    `two_point_correlator` for the first time, a dispatch that could not
    happen before (only n_uc>2 produced one). The polarized state fixes
    both exactly: <Sz>=1/2 on every site, <Sz Sz>=1/4 at every distance."""
    ic = _polarized_chain(n_uc, [(2, J2), (3, J3)])
    ic.gs_energy()
    assert getattr(ic._vumps_result, "multisite", False)
    assert ic.vev("Sz", 0).real == pytest.approx(0.5, abs=1e-9)
    for r in (1, 2, 3):
        assert ic.correlator("Sz", 0, "Sz", r).real == pytest.approx(0.25, abs=1e-9)


def test_finite_window_kpm_runs_on_a_long_range_chain():
    """`kpm_finite` builds its own finite window from the same canonicalized
    terms, so it has to survive a long-range one end to end -- not just
    `_window_hamiltonian` in isolation (checked separately below)."""
    ic = _polarized_chain(1, [(2, J2), (3, J3)])
    ic.gs_energy()
    out = np.array(ic.kpm_finite("Sz", 0, "Sz", 1, 6, delta=0.3))
    assert out.shape[0] == 2 and out.shape[1] > 0
    assert np.all(np.isfinite(out))


def test_reach_is_what_selects_the_sequential_solver():
    """A reach>1 Hamiltonian must go to the multi-site solver even at
    n_uc<=2, where the grouped path would otherwise take it."""
    ic = _polarized_chain(1, [(2, J2)])
    ic.gs_energy()
    assert getattr(ic._vumps_result, "multisite", False)
    # ... and a reach-1 one on the same cell must NOT.
    ic = _polarized_chain(1, [(1, J2)])
    ic.gs_energy()
    assert not getattr(ic._vumps_result, "multisite", False)


# -- 2. cell-size invariance -------------------------------------------------

def test_energy_is_independent_of_how_the_cell_is_sliced():
    """The same chain as a reach-2 coupling on a 1-site cell (the new
    sequential route) and as a reach-1 one on a 2-site cell (the
    already-validated grouped route): two code paths, one physical
    model."""
    ic1, ic2 = _ising_j2_chain(1), _ising_j2_chain(2)
    assert (ic1._reach_cells, ic2._reach_cells) == (2, 1)
    e1, e2 = ic1.gs_energy(), ic2.gs_energy()
    assert ic1.converged and ic2.converged
    assert getattr(ic1._vumps_result, "multisite", False)      # sequential
    assert not getattr(ic2._vumps_result, "multisite", False)  # grouped
    assert e1 == pytest.approx(e2, abs=1e-7)


@pytest.mark.parametrize("backend", BACKENDS)
def test_energy_is_independent_of_the_backend(backend):
    """The C++ backend has the same dispatch, so a reach-2 chain on a
    1-site cell must give the same energy density there as on the pure-
    Python one. `converged` is deliberately NOT asserted for
    itensor_version=3: its VUMPS Lanczos still stops on the eigenVALUE, so
    the gauge mismatch floors at ~1e-6 and the flag stays False even
    though the energy is right -- a known, separate C++ gap (see
    vumps.py's own "Convergence robustness" docstring section)."""
    e_ref = _ising_j2_chain(1).gs_energy()
    ic = _ising_j2_chain(1, backend=backend)
    assert ic._reach_cells == 2
    assert ic.gs_energy() == pytest.approx(e_ref, abs=1e-9)


@pytest.mark.skipif(3 not in BACKENDS, reason="mpscpp3 not compiled")
def test_v3_static_observables_say_why_they_cannot_answer():
    """itensor_version=3 answers the ENERGY for a long-range chain but not
    vev/correlator: `Chain::vumps_onsite_expectation` reads the grouped
    snapshot, and `vms_ground_state` -- which is what actually ran -- has
    never had a static-observable port. That is a pre-existing gap of the
    C++ sequential path (previously reachable only at n_uc>2), and the
    error has to name it rather than implying gs_energy was never
    called."""
    ic = _ising_j2_chain(1, backend=3)
    ic.gs_energy()
    with pytest.raises(RuntimeError, match="SEQUENTIAL multi-site solver"):
        ic.vev("Sz", 0)


# -- 3. canonicalization and gating ------------------------------------------

@pytest.mark.parametrize("n_uc", [1, 2, 3])
def test_canonicalization_translates_terms_without_losing_them(n_uc):
    """Every term keeps its coefficient, its operator names and its site
    SPACINGS -- canonicalization only translates, so a bond can neither be
    double-counted nor dropped."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"] * n_uc)
    h = (ic.get_operator("Sz", 0, group=-2) * ic.get_operator("Sz", 0, group=1)
         + ic.SzL[0] * ic.SzR[0]
         + ic.get_operator("Sz", 0, group=3) * ic.get_operator("Sz", 0, group=5)
         + 2.5 * ic.SzC[0])
    intra, inter = _canonicalize_hamiltonian(h, n_uc)
    got = []
    for term in list(intra.op) + list(inter.op):
        sites = [s for _n, s in term[1:]]
        assert min(sites) >= 0 and min(sites) < n_uc   # anchored on cell 0
        got.append((term[0], tuple(n for n, _s in term[1:]),
                    max(sites) - min(sites)))
    assert sorted(got) == sorted([
        (1.0, ("Sz", "Sz"), 3 * n_uc),      # group -2 -> 1
        (1.0, ("Sz", "Sz"), 2 * n_uc),      # L -> R
        (1.0, ("Sz", "Sz"), 2 * n_uc),      # group 3 -> 5
        (2.5, ("Sz",), 0),
    ])


def test_get_operator_accepts_integer_cell_offsets():
    ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"])
    assert ic.get_operator("Sz", 1, group=0).op == ic.get_operator("Sz", 1, "C").op
    assert ic.get_operator("Sz", 1, group=1).op == ic.get_operator("Sz", 1, "R").op
    assert ic.get_operator("Sz", 1, group=-1).op == ic.get_operator("Sz", 1, "L").op
    assert ic.get_operator("Sz", 0, group=3).op[0][1] == ["Sz", 6]
    with pytest.raises(ValueError):
        ic.get_operator("Sz", 0, group="X")


def test_finite_window_tiling_keeps_long_range_terms_inside_the_window():
    """`_window_hamiltonian` (kpm_finite's open-boundary approximation)
    tiles each inter-cell term only where its own copy still fits -- a
    reach-1 term drops one copy, a reach-2 term drops two."""
    n_uc, n_window = 1, 5
    intra, inter = _canonicalize_hamiltonian(
        infinitechain.Infinite_Spin_Chain(["1/2"]).get_operator("Sz", 0)
        * infinitechain.Infinite_Spin_Chain(["1/2"]).get_operator("Sz", 0, group=2),
        n_uc)
    window = _window_hamiltonian(intra, inter, n_uc, n_window)
    pairs = sorted(tuple(s for _n, s in term[1:]) for term in window.op)
    assert pairs == [(0, 2), (1, 3), (2, 4)]
    assert all(max(p) < n_window * n_uc for p in pairs)


def test_excitation_ansatz_rejects_a_long_range_hamiltonian():
    """The one algorithm that genuinely cannot do reach>1 says so, and
    says the ground state still works."""
    ic = _polarized_chain(1, [(2, J2)])
    ic.gs_energy()
    with pytest.raises(NotImplementedError, match="reaching 2 unit cells"):
        ic.excitation_energies(0.0)
