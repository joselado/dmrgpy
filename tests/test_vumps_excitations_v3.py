"""Coverage for the ITensor v3 C++ port of the tangent-space/quasiparticle
excitation ansatz (`mpscpp3/chain_session.h`'s
`Chain::vumps_excitation_energies`, wired into
`Infinite_Many_Body_Chain.excitation_energies`/`excitation_gap` via
`itensor_version=3`, `gs_method="vumps"`) -- cross-checked directly
against `itensor_version="python"`'s own `pyitensor/idmrg_excitations.py`
(already validated against an independently-converged MPSKit.jl D=2 TFIM
state, see that module's own "History" docstring section), since there is
no independent pytest coverage of the "python" excitation ansatz itself to
mirror here (only `test_vumps.py::test_gs_method_vumps_excitation_gap_works`,
a D=1 field-only smoke test).

Skipped automatically if mpscpp3 isn't compiled.
"""
import numpy as np
import pytest

from dmrgpy import cppext
from dmrgpy import infinitechain
from dmrgpy.pyitensor import vumps as pyvumps
from dmrgpy.pyitensor import idmrg_excitations as pyexc

pytestmark = pytest.mark.skipif(
    not cppext.available(3), reason="requires the compiled mpscpp3 (ITensor v3) extension")


def _tfim_chain(g, itensor_version):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=itensor_version)
    h = 4.0 * ic.SxC[0] * ic.SxR[0] + 2.0 * g * ic.SzC[0]
    ic.set_hamiltonian(h)
    return ic


def test_field_only_dispersion_is_flat():
    """D=1 field-polarized case (same as test_vumps.py's own
    test_gs_method_vumps_excitation_gap_works, run here through
    itensor_version=3 instead): E(k)=B for every k -- a single spin flip
    costs exactly the field energy, independent of momentum."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=3)
    ic.gs_method = "vumps"
    ic.maxm = 1
    ic.set_hamiltonian(2.0 * ic.SzC[0])
    for k in np.linspace(-np.pi, np.pi, 5):
        assert ic.excitation_energies(k, n=1)[0] == pytest.approx(2.0, abs=1e-6)
    assert ic.excitation_gap() == pytest.approx(2.0, abs=1e-6)


@pytest.mark.parametrize("D", [1, 2, 3])
def test_tfim_dispersion_matches_python_backend(D):
    """Direct cross-backend agreement across a momentum scan --
    Chain::vumps_excitation_energies is a line-for-line port of
    pyitensor/idmrg_excitations.py's own excitation_energies against the
    SAME converged VUMPSResult-equivalent mixed gauge, so the two must
    agree closely (TFIM is gapped away from g=1, so both the ground state
    and the excitation ansatz converge tightly and reproducibly at these
    D -- see test_vumps_v3.py's own analogous ground-state check). D=1 is
    a documented special case on BOTH backends (vumps.py's own D=1 mixed
    gauge never drives gauge_mismatch below a generic tol -- a phase-
    convention artifact of a 1x1 "isometry", not a bad fixed point, see
    test_vumps.py's own D=1 tests, which never assert .converged for the
    excitation-ansatz-adjacent cases either): `.converged` is therefore
    only asserted for D>1 here, though the dispersion values themselves
    still agree tightly regardless (checked below for every D)."""
    g = 1.5
    ic_v3 = _tfim_chain(g, 3)
    ic_v3.gs_method = "vumps"
    ic_v3.maxm = D
    ic_v3.vumps_nrestarts = 6
    ic_v3.gs_energy()
    if D > 1:
        assert ic_v3.converged

    ic_py = _tfim_chain(g, "python")
    result_py = pyvumps.vumps_ground_state(
        ic_py.site_types, ic_py._h_intra.op, ic_py._h_inter.op, ic_py.n_uc,
        D=D, tol=1e-10, maxiter=400, nrestarts=6)
    env_py = pyexc.build_excitation_environment(result_py)

    for k in np.linspace(-np.pi, np.pi, 9):
        e_v3 = ic_v3.excitation_energies(k, n=1)[0]
        e_py = pyexc.excitation_energies(env_py, k, n=1)[0]
        assert e_v3 == pytest.approx(e_py, abs=1e-6)


def test_excitation_gap_matches_python_backend():
    g = 1.5
    D = 2
    ic_v3 = _tfim_chain(g, 3)
    ic_v3.gs_method = "vumps"
    ic_v3.maxm = D
    ic_v3.vumps_nrestarts = 6

    ic_py = _tfim_chain(g, "python")
    ic_py.gs_method = "vumps"
    ic_py.maxm = D
    ic_py.vumps_nrestarts = 6

    gap_v3 = ic_v3.excitation_gap()
    gap_py = ic_py.excitation_gap()
    assert gap_v3 == pytest.approx(gap_py, abs=1e-4)


def test_n_uc2_heisenberg_dispersion_matches_python_backend():
    """n_uc=2 grouping path (dimerized unit cell) -- gapless/critical, so
    only a loose tolerance is asserted (both backends' non-convex restart
    searches can land on slightly different local optima -- see
    ROADMAP.md's own note on this), unlike the tight TFIM check above."""
    D = 2

    def make(itensor_version):
        ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"], itensor_version=itensor_version)
        h = (ic.SxC[0] * ic.SxC[1] + ic.SyC[0] * ic.SyC[1] + ic.SzC[0] * ic.SzC[1]
             + ic.SxC[1] * ic.SxR[0] + ic.SyC[1] * ic.SyR[0] + ic.SzC[1] * ic.SzR[0])
        ic.set_hamiltonian(h)
        ic.gs_method = "vumps"
        ic.maxm = D
        ic.vumps_nrestarts = 6
        return ic

    ic_v3 = make(3)
    ic_v3.gs_energy()
    ic_py = make("python")
    ic_py.gs_energy()

    for k in (0.0, np.pi / 2, np.pi):
        e_v3 = ic_v3.excitation_energies(k, n=1)[0]
        e_py = ic_py.excitation_energies(k, n=1)[0]
        assert e_v3 == pytest.approx(e_py, abs=5e-2)


def test_excitation_energies_requires_gs_method_vumps():
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=3)
    ic.set_hamiltonian(2.0 * ic.SzC[0])
    assert ic.gs_method == "idmrg"
    with pytest.raises(NotImplementedError):
        ic.excitation_energies(0.0)


def test_switching_gs_method_after_idmrg_run_still_works():
    """Regression check for the _session3_has_vumps bookkeeping: running
    gs_energy() once with the default gs_method="idmrg" (populating
    self._session3 with an idmrg-only snapshot), then switching to
    gs_method="vumps" and calling excitation_energies directly, must
    transparently rerun gs_energy() with the new gs_method rather than
    reusing the stale idmrg-only session (which has no VUMPS snapshot at
    all -- Chain::vumps_excitation_energies would otherwise raise
    "called before vumps_ground_state")."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=3)
    ic.set_hamiltonian(2.0 * ic.SzC[0])
    ic.maxm = 1
    ic.gs_energy()  # gs_method="idmrg" (default)
    ic.gs_method = "vumps"
    assert ic.excitation_energies(0.0, n=1)[0] == pytest.approx(2.0, abs=1e-6)
