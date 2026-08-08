"""Coverage for pyitensor/vumps.py's static correlator machinery
(onsite_expectation/two_point_correlator, computed directly from a
converged VUMPSResult's mixed-gauge {AC, AR} -- see vumps.py's own
"Static correlators" section for the derivation, following Vanderstraeten,
Haegeman, Verstraete, "Tangent-space methods for uniform matrix product
states", arXiv:1810.07006, Eq.(34)/(37)-(39)) and its
Infinite_Many_Body_Chain wiring (gs_method="vumps", vev/correlator).

Cross-checks used, in order of how rigorous they are: (1) exactly solvable
D=1 cases where the correlator is known in closed form -- a field-polarized
product state (every two-point correlator factorizes trivially, including
r=0 via Sz^2=1/4 for a spin-1/2 eigenstate) and a decoupled Heisenberg
singlet dimer (<Sz_i>=0, <Sz_0 Sz_1>=-1/4 exactly by SU(2) symmetry since
<S_0.S_1>=-3/4 for a singlet, and the two independent cells are exactly
UNcorrelated with each other, <Sz(cell m, site1) Sz(cell m+1, site0)>=0,
since the ground state literally factorizes across cells); (2) direct
agreement with pyitensor.idmrg's own independently-derived (dominant-
right-fixed-point-eigenproblem-based) onsite_expectation/
two_point_correlator on the SAME physical model at D>1 (TFIM) -- the two
backends' correlator machinery share no code beyond the low-level
idmrg._apply_transfer_from_left primitive, so agreement is a genuine
cross-check, not a shared-bug blind spot.
"""
import pytest

from dmrgpy import infinitechain
from dmrgpy.pyitensor import vumps


def _field_chain(B):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
    ic.gs_method = "vumps"
    ic.maxm = 1
    ic.set_hamiltonian(-B * ic.SzC[0])
    return ic


def _dimer_chain():
    """n_uc=2, intra-cell Heisenberg bond only (no inter-cell coupling at
    all) -- the exactly-solvable decoupled-singlet-dimer limit, D=1 after
    grouping, energy -3/4 per unit cell (see test_vumps.py's own
    _dimer_chain, which this mirrors)."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"], itensor_version="python")
    ic.gs_method = "vumps"
    ic.maxm = 1
    h = ic.SxC[0] * ic.SxC[1] + ic.SyC[0] * ic.SyC[1] + ic.SzC[0] * ic.SzC[1]
    ic.set_hamiltonian(h)
    return ic


def _tfim_chain(g):
    """H = -4*SxC[i]*SxC[i+1] - 2*g*SzC[i], same convention as
    test_vumps.py's own _tfim_chain -- the standard Pauli TFIM
    H=-sigma^x sigma^x - g*sigma^z."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
    h = -4.0 * ic.SxC[0] * ic.SxR[0] - 2.0 * g * ic.SzC[0]
    ic.set_hamiltonian(h)
    return ic


# == Exact D=1 cross-checks: field-polarized product state ==================

def test_field_polarized_onsite_matches_exact():
    ic = _field_chain(0.7)
    ic.gs_energy()
    assert ic.converged
    assert ic.vev("Sz", 0) == pytest.approx(0.5, abs=1e-8)


def test_field_polarized_onsite_negative_sign():
    ic = _field_chain(-0.4)
    ic.gs_energy()
    assert ic.vev("Sz", 0) == pytest.approx(-0.5, abs=1e-8)


@pytest.mark.parametrize("r", [0, 1, 2, 5])
def test_field_polarized_correlator_factorizes(r):
    """A product state has no correlations at all: <Sz(0)Sz(r)> =
    <Sz>^2 = 0.25 for every r>0, and Sz^2=1/4 exactly at r=0 too (a
    spin-1/2 eigenstate) -- so the SAME value 0.25 holds for every r,
    including across several unit-cell boundaries (n_uc=1 here, so every
    r>=1 already crosses at least one)."""
    ic = _field_chain(0.7)
    ic.gs_energy()
    assert ic.correlator("Sz", 0, "Sz", r) == pytest.approx(0.25, abs=1e-8)


# == Exact D=1 cross-checks: decoupled singlet dimer =========================

def test_dimer_onsite_is_zero():
    """A singlet has zero magnetization on both sites."""
    ic = _dimer_chain()
    ic.gs_energy()
    assert ic.converged
    assert ic.vev("Sz", 0) == pytest.approx(0.0, abs=1e-8)
    assert ic.vev("Sz", 1) == pytest.approx(0.0, abs=1e-8)


def test_dimer_intracell_correlator_matches_singlet_value():
    """<S_0.S_1> = -3/4 for an exact singlet; by SU(2) symmetry each
    Cartesian component contributes equally, so <Sz_0 Sz_1> = -1/4 --
    both operators land in the SAME AC (grouped-supersite) cell here
    (cell_offset=0 in two_point_correlator's own internal split)."""
    ic = _dimer_chain()
    ic.gs_energy()
    assert ic.correlator("Sz", 0, "Sz", 1) == pytest.approx(-0.25, abs=1e-6)


def test_dimer_intercell_correlator_is_exactly_zero():
    """Site 1 of one cell and site 0 of the next cell are in DIFFERENT,
    completely decoupled singlets -- exactly uncorrelated. This is the
    ONLY case among these dimer tests that exercises
    two_point_correlator's cross-cell branch (cell_offset=1, propagating
    through the AR transfer tensor rather than staying inside one AC
    cell)."""
    ic = _dimer_chain()
    ic.gs_energy()
    assert ic.correlator("Sz", 1, "Sz", 1) == pytest.approx(0.0, abs=1e-6)


def test_dimer_same_site_correlator_is_spin_squared():
    """<Sz_i^2> = 1/4 for a spin-1/2 eigenbasis, regardless of state --
    exercises two_point_correlator's r=0 (M_j @ M_i) same-site branch."""
    ic = _dimer_chain()
    ic.gs_energy()
    assert ic.correlator("Sz", 0, "Sz", 0) == pytest.approx(0.25, abs=1e-6)


def test_dimer_correlator_across_two_full_cells_is_zero():
    """r=3 from p_i=0 lands at cell_offset=1 (site 3 = cell 1, site 1) --
    still a different, decoupled singlet from cell 0's own site 0, so
    still exactly uncorrelated. Distinct from the r=1 cross-cell case
    above (different p_i/p_j and a different absolute r), so this is not
    a duplicate of test_dimer_intercell_correlator_is_exactly_zero."""
    ic = _dimer_chain()
    ic.gs_energy()
    assert ic.correlator("Sz", 0, "Sz", 3) == pytest.approx(0.0, abs=1e-6)


# == Direct vumps.py module-level functions (no Infinite_Many_Body_Chain) ===

def test_vumps_onsite_expectation_function_matches_exact():
    ic = _field_chain(1.0)
    result = vumps.vumps_ground_state(ic.site_types, ic._h_intra.op, ic._h_inter.op,
                                       ic.n_uc, D=1, maxiter=50)
    assert vumps.onsite_expectation(result, "Sz", 0) == pytest.approx(0.5, abs=1e-8)


def test_vumps_two_point_correlator_function_matches_exact():
    ic = _dimer_chain()
    result = vumps.vumps_ground_state(ic.site_types, ic._h_intra.op, ic._h_inter.op,
                                       ic.n_uc, D=1, maxiter=50)
    assert vumps.two_point_correlator(result, "Sz", 0, "Sz", 1) == pytest.approx(-0.25, abs=1e-6)


def test_vumps_onsite_expectation_rejects_p_out_of_range():
    ic = _field_chain(1.0)
    result = vumps.vumps_ground_state(ic.site_types, ic._h_intra.op, ic._h_inter.op,
                                       ic.n_uc, D=1, maxiter=50)
    with pytest.raises(ValueError):
        vumps.onsite_expectation(result, "Sz", 1)


def test_vumps_two_point_correlator_rejects_negative_r():
    ic = _field_chain(1.0)
    result = vumps.vumps_ground_state(ic.site_types, ic._h_intra.op, ic._h_inter.op,
                                       ic.n_uc, D=1, maxiter=50)
    with pytest.raises(ValueError):
        vumps.two_point_correlator(result, "Sz", 0, "Sz", -1)


# == Cross-check against idmrg.py's independently-derived correlator machinery ==

@pytest.mark.parametrize("D", [2, 3])
def test_vumps_correlator_matches_idmrg_backend_on_tfim(D):
    """Same physical model (TFIM, gapped, away from criticality), two
    independently-derived correlator implementations (idmrg.py's
    dominant-right-fixed-point eigenproblem vs vumps.py's exact AC/AR
    mixed-gauge contraction) converged independently at the same target
    bond dimension -- agreement confirms neither has a sign/convention
    bug the other doesn't share. Looser tolerance than the exact D=1
    cases above since both are independent NUMERICAL optimizations
    (growing-algorithm truncation vs VUMPS's own restart search) that
    need not land on bit-identical states even at the same D. `.converged`
    itself is NOT asserted here at D=3 -- mirrors
    tests/test_vumps_v3.py's own test_tfim_matches_python_backend, whose
    docstring documents that VUMPS's gauge_mismatch convergence tail at
    D=3 can stay above `tol` for many more outer iterations than the
    ENERGY itself needs to settle (confirmed there: 5 independent D=3
    runs at maxiter=800 all landed on the identical energy to 9 digits
    while `.converged` stayed False) -- so the correlator VALUES, not the
    convergence flag, are the thing actually being checked."""
    g = 1.5

    ic_v = _tfim_chain(g)
    ic_v.gs_method = "vumps"
    ic_v.maxm = D
    ic_v.maxiter = 800
    ic_v.vumps_nrestarts = 6
    ic_v.gs_energy()

    ic_i = _tfim_chain(g)
    ic_i.gs_method = "idmrg"
    ic_i.maxm = D
    ic_i.maxiter = 200
    ic_i.gs_energy()

    sz_v = ic_v.vev("Sz", 0)
    sz_i = ic_i.vev("Sz", 0)
    assert sz_v == pytest.approx(sz_i, abs=0.05)

    for r in (1, 2, 3):
        corr_v = ic_v.correlator("Sz", 0, "Sz", r)
        corr_i = ic_i.correlator("Sz", 0, "Sz", r)
        assert corr_v == pytest.approx(corr_i, abs=0.05)
