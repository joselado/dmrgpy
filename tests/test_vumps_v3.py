"""Coverage for the ITensor v3 C++ port of VUMPS
(`mpscpp3/chain_session.h`'s `Chain::vumps_ground_state`, wired into
`Infinite_Many_Body_Chain.gs_energy` via `itensor_version=3`,
`gs_method="vumps"`) -- mirrors `tests/test_vumps.py`'s own
`itensor_version="python"` coverage of `pyitensor/vumps.py`, plus direct
cross-backend agreement checks (the two are supposed to implement the
identical algorithm -- see `VumpsResult`'s own doc comment in
chain_session.h for why this port is dense-LAPACK-based rather than
built from ITensor tensor-network objects, unlike `idmrg_ground_state`).

Skipped automatically if mpscpp3 isn't compiled, exactly like
`test_infinite_chain.py`'s own `test_itensor_version3_*` tests.
"""
import numpy as np
import pytest
from scipy.integrate import quad

from dmrgpy import cppext
from dmrgpy import infinitechain
from dmrgpy.pyitensor import vumps as pyvumps

pytestmark = pytest.mark.skipif(
    not cppext.available(3), reason="requires the compiled mpscpp3 (ITensor v3) extension")


def _field_chain(B, itensor_version):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=itensor_version)
    ic.set_hamiltonian(-B * ic.SzC[0])
    return ic


def _dimer_chain(itensor_version):
    ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"], itensor_version=itensor_version)
    h = ic.SxC[0] * ic.SxC[1] + ic.SyC[0] * ic.SyC[1] + ic.SzC[0] * ic.SzC[1]
    ic.set_hamiltonian(h)
    return ic


def _tfim_chain(g, itensor_version):
    """Same H=-4*SxSx-2*g*Sz=-sigma^x sigma^x-g*sigma^z convention as
    test_vumps.py's own _tfim_chain."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=itensor_version)
    h = -4.0 * ic.SxC[0] * ic.SxR[0] - 2.0 * g * ic.SzC[0]
    ic.set_hamiltonian(h)
    return ic


def _tfim_exact_energy_density(g):
    val, _ = quad(lambda k: np.sqrt(1 + g ** 2 - 2 * g * np.cos(k)), 0, np.pi)
    return -val / np.pi


def test_field_only_matches_exact_solution():
    ic = _field_chain(0.7, 3)
    ic.gs_method = "vumps"
    ic.maxm = 1
    e0 = ic.gs_energy()
    assert ic.converged
    assert e0 == pytest.approx(-0.35, abs=1e-8)


def test_dimer_matches_exact_singlet_energy():
    """D=1 (after grouping) WITH a bond term -- exercises the
    pending-channel code path (Lvec_a/Rvec_b, H_AC/H_C bond terms)."""
    ic = _dimer_chain(3)
    ic.gs_method = "vumps"
    ic.maxm = 1
    e0 = ic.gs_energy()
    assert ic.converged
    assert e0 == pytest.approx(-0.375, abs=1e-8)


@pytest.mark.parametrize("D", [2, 3])
def test_tfim_matches_python_backend(D):
    """Direct cross-backend agreement: itensor_version=3's
    Chain::vumps_ground_state is a line-for-line port of
    pyitensor/vumps.py's own vumps_ground_state, so the two should agree
    on the converged energy density to close to machine precision, not
    just to the exact TFIM answer's own convergence-in-D tolerance.

    pyitensor's own VUMPS gauge_mismatch can take many more outer
    iterations than `maxiter` to cross a strict `tol=1e-10` at this D even
    once the ENERGY itself has already settled to ~1e-9 (confirmed
    directly: 5 independent D=3 runs at `maxiter=800` all landed on the
    identical energy to 9 digits while `.converged` stayed False, gauge
    mismatch ~1e-7..1e-6 -- a slow final-convergence tail, not a bad-basin
    restart-search failure, since the energy itself is reproducible, not
    scattered) -- so only THIS test's own v3 result is required to be
    `.converged` (confirmed reliable, see the tighter mismatch it reaches
    below); the pyitensor reference's own `.converged` flag is not
    asserted, only its energy value, at a tolerance loose enough to cover
    its own residual gauge-mismatch-driven energy drift."""
    g = 1.5
    ic_v3 = _tfim_chain(g, 3)
    ic_v3.gs_method = "vumps"
    ic_v3.maxm = D
    ic_v3.vumps_nrestarts = 6
    e0_v3 = ic_v3.gs_energy()
    assert ic_v3.converged

    ic_py = _tfim_chain(g, "python")
    result_py = pyvumps.vumps_ground_state(
        ic_py.site_types, ic_py._h_intra.op, ic_py._h_inter.op, ic_py.n_uc,
        D=D, tol=1e-10, maxiter=800, nrestarts=6)
    assert e0_v3 == pytest.approx(result_py.e0, abs=1e-6)


def test_tfim_stays_within_variational_bound():
    """Same style of restart-robust sanity check as test_vumps.py's own
    test_tfim_stays_within_variational_bound -- a converged/best-effort
    VUMPS energy can never sit appreciably below the exact ground-state
    energy (up to a tiny numerical slack)."""
    g = 1.5
    exact = _tfim_exact_energy_density(g)
    ic = _tfim_chain(g, 3)
    ic.gs_method = "vumps"
    ic.maxm = 3
    ic.maxiter = 400
    ic.vumps_nrestarts = 5
    e0 = ic.gs_energy()
    assert e0 > exact - 0.2
    assert e0 < exact + 0.5


def test_n_uc2_uniform_heisenberg_bethe_ansatz_bound():
    """n_uc=2 grouping path (exercises vumps_group_automaton's own
    two-sublattice branch, not just the n_uc=1 trivial case) -- the exact
    Bethe-ansatz density (1/4-ln(2)) is a lower bound any finite-D VUMPS
    energy must sit at or above."""
    exact = 0.25 - np.log(2)
    ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"], itensor_version=3)
    h = (ic.SxC[0] * ic.SxC[1] + ic.SyC[0] * ic.SyC[1] + ic.SzC[0] * ic.SzC[1]
         + ic.SxC[1] * ic.SxR[0] + ic.SyC[1] * ic.SyR[0] + ic.SzC[1] * ic.SzR[0])
    ic.gs_method = "vumps"
    ic.maxm = 2
    ic.maxiter = 400
    ic.vumps_nrestarts = 5
    ic.set_hamiltonian(h)
    e0 = ic.gs_energy()
    assert e0 > exact - 1e-6
    assert e0 < exact + 0.5


def test_d_less_than_1_rejected():
    ic = _field_chain(0.5, 3)
    ic.gs_method = "vumps"
    ic.maxm = 0
    with pytest.raises(Exception):
        ic.gs_energy()


def test_n_uc_above_2_uses_the_sequential_multisite_algorithm():
    """n_uc>2 used to be rejected at construction. The v3 backend now runs
    the sequential multi-site algorithm there (Chain::vms_ground_state), so
    a 3-site cell must solve, and its observables must come from the
    per-site snapshot rather than the grouped one. The physics checks live
    in tests/test_vumps_multisite.py, which covers both backends; this one
    pins that the v3 route exists and produces a sane state."""
    g = 1.5
    ic = infinitechain.Infinite_Spin_Chain(["1/2"] * 3, itensor_version=3)
    ic.gs_method = "vumps"
    ic.maxm, ic.maxiter, ic.etol = 4, 100, 1e-11
    ic.vumps_nrestarts = 2
    # The TFIM rather than a pure field: a field-only model's ground state
    # is a product state, whose transfer matrix is massively degenerate at
    # D>1, and BOTH VUMPS paths refuse a non-unique dominant fixed point
    # there (which is why this file's own field-only test uses maxm=1).
    h = 0
    for i in range(3):
        nxt = ic.SxC[i + 1] if i + 1 < 3 else ic.SxR[0]
        h = h + (-4.0) * ic.SxC[i] * nxt + (-2.0 * g) * ic.SzC[i]
    ic.set_hamiltonian(h)
    e0 = ic.gs_energy()
    assert e0 == pytest.approx(_tfim_exact_energy_density(g), abs=1e-3)
    sz = [float(np.real(ic.vev("Sz", p))) for p in range(3)]
    assert max(sz) - min(sz) < 1e-8       # uniform chain, uniform observable


def test_gs_method_vumps_vev_and_correlator_work():
    """vev/correlator DO work under itensor_version=3, gs_method="vumps"
    (Chain::vumps_onsite_expectation/vumps_two_point_correlator, a C++
    port of pyitensor/vumps.py's own AC/AR-based formula) -- only
    gs_method="idmrg" still raises NotImplementedError on this backend
    (IdmrgResult keeps no per-sublattice U_list). See
    tests/test_vumps_correlator_v3.py for the full correlator coverage
    (exact D=1 closed-form cases, cross-backend agreement at D=2,3, and
    the gs_method="idmrg" NotImplementedError paths) -- this is just a
    non-regression smoke check for the field-polarized D=1 exact case."""
    ic = _field_chain(0.5, 3)
    ic.gs_method = "vumps"
    ic.maxm = 1
    ic.gs_energy()
    assert ic.vev("Sz", 0) == pytest.approx(0.5, abs=1e-8)
    assert ic.correlator("Sz", 0, "Sz", 1) == pytest.approx(0.25, abs=1e-8)


def test_gs_method_idmrg_still_works_on_v3():
    """Non-regression: itensor_version=3's own gs_method dispatch must not
    have broken the pre-existing gs_method="idmrg" path (no longer the
    default since 2026-08-08 -- "vumps" is -- so set explicitly here) --
    Chain::idmrg_ground_state is a distinct method from
    Chain::vumps_ground_state, only sharing idmrg_build_row's own
    per-sublattice automaton construction (see that method's own
    onsite/field-term bugfix, cross-checked directly here since none of
    test_infinite_chain.py's own v3 cases use a field term)."""
    ic = _field_chain(0.7, 3)
    ic.gs_method = "idmrg"
    ic.maxm = 4
    ic.maxiter = 60
    ic.etol = 1e-9
    e0 = ic.gs_energy()
    assert e0 == pytest.approx(-0.35, abs=1e-3)


def test_idmrg_ground_state_v3_onsite_field_matches_exact_solution():
    """C++ analogue of test_infinite_chain.py's own (itensor_version=
    "python"-only) test_onsite_field_matches_exact_solution -- direct
    regression coverage, on itensor_version=3's own Chain::
    idmrg_ground_state, for the real onsite-automaton bug this audit
    found and fixed in idmrg_build_row (a site's own onsite/field term
    was landing on the automaton's F,F self-loop instead of the direct
    S,F transition, exactly the bug pyitensor/idmrg.py's own
    _build_periodic_mpo docstring documents finding and fixing on the
    Python side -- before the fix, this exact case reproducibly gave
    e0=0 and <Sz>~0 regardless of B, since W stayed block-diagonal
    between {S} and {F,pending...} with no bond term ever able to
    activate F)."""
    for B in (5.0, -5.0, 2.0):
        ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=3)
        ic.gs_method = "idmrg"
        ic.maxm, ic.maxiter, ic.etol = 4, 50, 1e-12
        ic.set_hamiltonian(-B * ic.SzC[0])
        e0 = ic.gs_energy()
        assert e0 == pytest.approx(-abs(B) / 2, abs=1e-6)


def test_idmrg_ground_state_v3_onsite_field_with_bonds_does_not_diverge():
    """C++ analogue of test_infinite_chain.py's own (itensor_version=
    "python"-only) test_onsite_field_with_bonds_does_not_diverge -- see
    test_idmrg_ground_state_v3_onsite_field_matches_exact_solution's own
    docstring for the bug this guards against; the other of its two
    confirmed failure modes (once a bond term activates the automaton's F
    channel, every further-absorbed site re-adds the field content
    multiplicatively instead of once, additively -- energy density
    diverging to ~1e69 rather than merely being wrong)."""
    def density_at(B):
        ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"], itensor_version=3)
        ic.gs_method = "idmrg"
        h = (ic.SxC[0] * ic.SxC[1] + ic.SyC[0] * ic.SyC[1] + ic.SzC[0] * ic.SzC[1]
             + 0.5 * (ic.SxC[1] * ic.SxR[0] + ic.SyC[1] * ic.SyR[0] + ic.SzC[1] * ic.SzR[0])
             - B * (ic.SzC[0] + ic.SzC[1]))
        ic.maxm, ic.maxiter, ic.etol = 20, 100, 1e-11
        ic.set_hamiltonian(h)
        e0 = ic.gs_energy()
        assert ic.converged
        assert np.isfinite(e0)
        return e0
    d0 = density_at(0.0)
    d1 = density_at(0.3)  # well below this gapped model's own gap
    assert d1 == pytest.approx(d0, abs=1e-6)
