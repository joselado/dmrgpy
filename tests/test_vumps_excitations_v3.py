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
    ic.gs_method = "idmrg"  # "vumps" (the default since 2026-08-08) DOES work -- see above
    with pytest.raises(NotImplementedError):
        ic.excitation_energies(0.0)


def test_switching_gs_method_after_idmrg_run_still_works():
    """Regression check for the _session3_has_vumps bookkeeping: running
    gs_energy() once with gs_method="idmrg" explicitly (populating
    self._session3 with an idmrg-only snapshot), then switching to
    gs_method="vumps" and calling excitation_energies directly, must
    transparently rerun gs_energy() with the new gs_method rather than
    reusing the stale idmrg-only session (which has no VUMPS snapshot at
    all -- Chain::vumps_excitation_energies would otherwise raise
    "called before vumps_ground_state")."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=3)
    ic.set_hamiltonian(2.0 * ic.SzC[0])
    ic.maxm = 1
    ic.gs_method = "idmrg"
    ic.gs_energy()
    ic.gs_method = "vumps"
    assert ic.excitation_energies(0.0, n=1)[0] == pytest.approx(2.0, abs=1e-6)


def _heisenberg_chain(D, itensor_version=3):
    ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"], itensor_version=itensor_version)
    h = (ic.SxC[0] * ic.SxC[1] + ic.SyC[0] * ic.SyC[1] + ic.SzC[0] * ic.SzC[1]
         + ic.SxC[1] * ic.SxR[0] + ic.SyC[1] * ic.SyR[0] + ic.SzC[1] * ic.SzR[0])
    ic.set_hamiltonian(h)
    ic.gs_method = "vumps"
    ic.maxm = D
    ic.vumps_nrestarts = 6
    return ic


@pytest.mark.parametrize("model,D,dim", [("tfim", 2, 4), ("tfim", 3, 9),
                                          ("tfim", 4, 16), ("heis", 2, 12),
                                          ("heis", 3, 27), ("heis", 4, 48)])
@pytest.mark.parametrize("nev", [1, 2, 3])
def test_lanczos_h_eff_matches_the_dense_solver(model, D, dim, nev):
    """Above `Chain::vumps_h_eff_dense_max_` the excitation ansatz solves
    H_eff(k) by Lanczos on its action instead of assembling it -- the two
    solvers must return the same energies, and the way to check that is on
    ONE converged state, with the solver forced, rather than across two
    runs whose ground states already differ.

    `dense_max` is what forces it (a static constexpr cannot be
    monkeypatched the way pyitensor's own `_DENSE_EIG_MAX` is): -1 keeps
    the built-in threshold, so every chain here (dim<=48) takes the dense
    path, and 0 puts the same chain on the Lanczos path.

    The n_uc=2 Heisenberg rows are the ones that matter, and they are here
    because they caught a real bug rather than to be thorough: their
    H_eff(k) spectrum is pairwise degenerate away from k=0, and a plain
    single-vector Lanczos finds only ONE copy of a degenerate eigenvalue,
    so asking it for the lowest three returned three DISTINCT eigenvalues
    (0.298598886, 0.299768235, 1.337377610 at k=0.37, D=2) where the dense
    answer is 0.298598886 twice and then 0.299768235 -- everything after
    the first copy shifted up, an error of 1.0 at nev=3, and every value
    returned a genuine eigenvalue, so no residual test sees it. What fixes
    it is the deflation and the per-run generic start vector in
    `vx_lanczos_lowest`, which is what these rows pin."""
    if nev >= dim:
        pytest.skip("nev must be below the dimension for the iterative path")
    ic = _tfim_chain(1.5, 3) if model == "tfim" else _heisenberg_chain(D)
    if model == "tfim":
        ic.gs_method = "vumps"
        ic.maxm = D
        ic.vumps_nrestarts = 6
    ic.gs_energy()
    session = ic._session3
    for k in list(np.linspace(-np.pi, np.pi, 5)) + [0.37]:
        dense = np.array(session.vumps_excitation_energies(k, nev, -1))
        lanczos = np.array(session.vumps_excitation_energies(k, nev, 0))
        assert len(lanczos) == len(dense)
        assert lanczos == pytest.approx(dense, abs=1e-8)


def test_momentum_resolvent_cache_does_not_leak_between_states():
    """The two channel resolvents of the momentum most recently asked for
    are cached on the Chain (`Chain::vumps_exc_resolvents`), so a chain
    that solves a new ground state and is then asked about the SAME
    momentum again must answer from the new state, not from the factors
    the old one left behind. A stale cache would show up here as the D=1
    field-only answer surviving into the D=2 chain."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=3)
    ic.gs_method = "vumps"
    ic.maxm = 1
    ic.set_hamiltonian(2.0 * ic.SzC[0])
    assert ic.excitation_energies(0.7, n=1)[0] == pytest.approx(2.0, abs=1e-6)

    ic2 = _tfim_chain(1.5, 3)
    ic2.gs_method = "vumps"
    ic2.maxm = 2
    ic2.vumps_nrestarts = 6
    e_first = ic2.excitation_energies(0.7, n=1)[0]
    # A second gs_energy() re-solves (infinitechain.py builds a fresh
    # Chain each time rather than caching), so this asks the same
    # momentum of a session that has never been asked it before, while
    # the momentum walk below asks one session for several momenta in two
    # different orders -- that walk is the direct test of the cache, this
    # one only says the answer does not depend on the session's history.
    ic2.gs_energy()
    assert ic2.excitation_energies(0.7, n=1)[0] == pytest.approx(e_first, abs=1e-6)
    # And the momentum scan itself must not depend on the order it is
    # walked in, which a resolvent held for the wrong momentum would break.
    ks = [0.0, 0.7, 1.9, np.pi]
    forward = [ic2.excitation_energies(k, n=1)[0] for k in ks]
    backward = [ic2.excitation_energies(k, n=1)[0] for k in reversed(ks)]
    assert forward == pytest.approx(list(reversed(backward)), abs=1e-10)
