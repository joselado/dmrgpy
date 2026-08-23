"""Coverage for the ITensor v3 C++ port of VUMPS's static-correlator
machinery (`mpscpp3/chain_session.h`'s `Chain::vumps_onsite_expectation`/
`Chain::vumps_two_point_correlator`, wired into
`Infinite_Many_Body_Chain.vev`/`correlator` via `itensor_version=3`,
`gs_method="vumps"`) -- mirrors `tests/test_vumps_correlator.py`'s own
`itensor_version="python"` coverage of `pyitensor/vumps.py`'s
`onsite_expectation`/`two_point_correlator` (same exact D=1 closed-form
cases: a field-polarized product state and a decoupled Heisenberg singlet
dimer), plus direct cross-backend agreement against `itensor_version=
"python"` at D>1 on TFIM (see `tests/test_vumps_excitations_v3.py`'s own
docstring for why this kind of direct comparison is a genuine cross-check,
not a shared-bug blind spot: the C++ port is a from-scratch dense-array
translation, sharing no code with pyitensor/vumps.py beyond the algorithm
itself).

Skipped automatically if mpscpp3 isn't compiled, exactly like
`test_vumps_v3.py`'s own tests.
"""
import pytest

from dmrgpy import cppext
from dmrgpy import infinitechain

pytestmark = pytest.mark.skipif(
    not cppext.available(3), reason="requires the compiled mpscpp3 (ITensor v3) extension")


def _field_chain(B, itensor_version=3):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=itensor_version)
    ic.gs_method = "vumps"
    ic.maxm = 1
    ic.set_hamiltonian(-B * ic.SzC[0])
    return ic


def _dimer_chain(itensor_version=3):
    """n_uc=2, intra-cell Heisenberg bond only (no inter-cell coupling at
    all) -- the exactly-solvable decoupled-singlet-dimer limit, D=1 after
    grouping, energy -3/4 per unit cell -- see test_vumps_correlator.py's
    own _dimer_chain, which this mirrors."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"], itensor_version=itensor_version)
    ic.gs_method = "vumps"
    ic.maxm = 1
    h = ic.SxC[0] * ic.SxC[1] + ic.SyC[0] * ic.SyC[1] + ic.SzC[0] * ic.SzC[1]
    ic.set_hamiltonian(h)
    return ic


def _tfim_chain(g, itensor_version):
    """H = -4*SxC[i]*SxC[i+1] - 2*g*SzC[i], same convention as
    test_vumps_v3.py's/test_vumps_correlator.py's own _tfim_chain."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=itensor_version)
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
    spin-1/2 eigenstate)."""
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
    (cell_offset=0 in Chain::vumps_two_point_correlator's own internal
    split)."""
    ic = _dimer_chain()
    ic.gs_energy()
    assert ic.correlator("Sz", 0, "Sz", 1) == pytest.approx(-0.25, abs=1e-6)


def test_dimer_intercell_correlator_is_exactly_zero():
    """Site 1 of one cell and site 0 of the next cell are in DIFFERENT,
    completely decoupled singlets -- exactly uncorrelated. The ONLY case
    among these dimer tests that exercises the cross-cell branch
    (cell_offset=1, propagating through the AR transfer tensor)."""
    ic = _dimer_chain()
    ic.gs_energy()
    assert ic.correlator("Sz", 1, "Sz", 1) == pytest.approx(0.0, abs=1e-6)


def test_dimer_same_site_correlator_is_spin_squared():
    """<Sz_i^2> = 1/4 for a spin-1/2 eigenbasis, regardless of state --
    exercises the r=0 (Mj@Mi) same-site branch."""
    ic = _dimer_chain()
    ic.gs_energy()
    assert ic.correlator("Sz", 0, "Sz", 0) == pytest.approx(0.25, abs=1e-6)


def test_dimer_correlator_across_two_full_cells_is_zero():
    """r=3 from p_i=0 lands at cell_offset=1 (site 3 = cell 1, site 1) --
    still a different, decoupled singlet from cell 0's own site 0."""
    ic = _dimer_chain()
    ic.gs_energy()
    assert ic.correlator("Sz", 0, "Sz", 3) == pytest.approx(0.0, abs=1e-6)


# == Direct cross-check against itensor_version="python" at D>1 (TFIM) ======

@pytest.mark.parametrize("D", [2, 3])
def test_correlator_matches_python_backend_on_tfim(D):
    """Same physical model, same algorithm (VUMPS), two independent
    implementations (this C++ port vs pyitensor/vumps.py) -- see
    test_vumps_v3.py's own test_tfim_matches_python_backend for why
    `.converged` itself is not asserted at D=3 (the gauge-mismatch
    convergence tail can lag well behind the energy/observables actually
    settling); the correlator VALUES are the thing checked here."""
    g = 1.5

    ic_v3 = _tfim_chain(g, 3)
    ic_v3.gs_method = "vumps"
    ic_v3.maxm = D
    ic_v3.maxiter = 800
    ic_v3.vumps_nrestarts = 6
    ic_v3.gs_energy()

    ic_py = _tfim_chain(g, "python")
    ic_py.gs_method = "vumps"
    ic_py.maxm = D
    ic_py.maxiter = 800
    ic_py.vumps_nrestarts = 6
    ic_py.gs_energy()

    assert ic_v3.e0 == pytest.approx(ic_py.e0, abs=1e-6)

    sz_v3 = ic_v3.vev("Sz", 0)
    sz_py = ic_py.vev("Sz", 0)
    assert sz_v3 == pytest.approx(sz_py, abs=1e-6)

    for r in (1, 2, 3):
        corr_v3 = ic_v3.correlator("Sz", 0, "Sz", r)
        corr_py = ic_py.correlator("Sz", 0, "Sz", r)
        assert corr_v3 == pytest.approx(corr_py, abs=1e-6)


# == Error paths ==============================================================

def test_vev_rejects_unknown_gs_method_on_v3():
    """Both "vumps" (this file) and "idmrg"
    (tests/test_idmrg_correlator_v3.py) support vev/correlator on the v3
    backend now; anything else must still raise."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=3)
    ic.set_hamiltonian(2.0 * ic.SzC[0])
    ic.gs_method = "nonsense"
    with pytest.raises(NotImplementedError):
        ic.vev("Sz", 0)


def test_correlator_rejects_unknown_gs_method_on_v3():
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=3)
    ic.set_hamiltonian(2.0 * ic.SzC[0])
    ic.gs_method = "nonsense"
    with pytest.raises(NotImplementedError):
        ic.correlator("Sz", 0, "Sz", 1)


def test_vev_rejects_p_out_of_range():
    ic = _field_chain(1.0)
    ic.gs_energy()
    with pytest.raises(ValueError):
        ic.vev("Sz", 1)


def test_correlator_rejects_negative_r():
    """r<0 is only checked Python-side by pyitensor/vumps.py's own
    two_point_correlator (ValueError); on itensor_version=3 it is instead
    caught by Chain::vumps_two_point_correlator's own C++-side ITError,
    which pybind11 surfaces as a plain RuntimeError -- mirrors
    test_vumps_v3.py's own convention of asserting `Exception` (not a
    specific subclass) for this backend's input-validation paths."""
    ic = _field_chain(1.0)
    ic.gs_energy()
    with pytest.raises(Exception):
        ic.correlator("Sz", 0, "Sz", -1)


def test_switching_gs_method_after_idmrg_run_still_works():
    """Regression check for the _session3_has_vumps bookkeeping (mirrors
    test_vumps_excitations_v3.py's own analogous test): running gs_energy()
    once with gs_method="idmrg" explicitly (populating self._session3 with
    an idmrg-only snapshot), then switching to gs_method="vumps" and
    calling vev/correlator directly, must transparently rerun gs_energy()
    with the new gs_method rather than reusing the stale idmrg-only
    session."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=3)
    ic.set_hamiltonian(-0.7 * ic.SzC[0])
    ic.maxm = 1
    ic.gs_method = "idmrg"
    ic.gs_energy()
    ic.gs_method = "vumps"
    assert ic.vev("Sz", 0) == pytest.approx(0.5, abs=1e-8)
