"""dmrg_method="CBE" (one-site DMRG with controlled bond expansion) on
itensor_version="julia_live", mpsjulialive/cbe.jl.

What is pinned here:

- the julia_live route runs cbe.jl's get_gs_cbe, not the two-site dmrg();
- CBE reaches ED on a Heisenberg chain (two random starts) and on a
  spinless fermion chain (Jordan-Wigner strings in the MPO), and agrees
  with itensor_version="python"'s CBE;
- from the Neel product state, all bonds 1, the expansion lets the bonds
  grow and the solve reaches the ground state;
- tevol_method="TDVP_CBE" (one-site TDVP with the same expansion) follows
  ED after a quench from the Neel state, takes complex time steps (the
  TDZ route), and gives the quench correlator of submode="TD".

Every case pins its sweep schedule; the Julia random starts are seeded
through Julia's own generator.
"""

import contextlib
import io
import types
import warnings

import numpy as np
import pytest

from dmrgpy import fermionchain, groundstate, spinchain, timedependent

from _helpers import julia_available

JULIA = [pytest.param("julia_live", id="julia_live")]


def _require_julia(version):
    ok, reason = julia_available()
    if not ok:
        pytest.skip("requires a working juliacall/Julia toolchain: %s" % reason)


def _quiet(f):
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return f()


def _seed(k):
    from dmrgpy.mpsjulialive.juliasession import Main
    Main.seval("import Random; Random.seed!(%d)" % k)


def _links(wf):
    from dmrgpy.mpsjulialive.juliasession import Main
    f = Main.seval("psi -> [dim(linkind(psi,i)) for i in 1:length(psi)-1]")
    return list(f(wf.jlmps))


def _heisenberg(n, version, **kw):
    sc = spinchain.Spin_Chain([2] * n, itensor_version=version, **kw)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]
    h = h + 0.3 * sc.Sz[0]
    sc.set_hamiltonian(h)
    return sc


class _Recorder:
    """Stands in for the Julia Main module inside
    mpsjulialive/groundstate.py and records which entry point it hands
    out."""

    def __init__(self, main):
        self._main, self.calls = main, []

    def __getattr__(self, name):
        self.calls.append(name)
        return getattr(self._main, name)


@pytest.mark.parametrize("version", JULIA)
def test_julia_cbe_heisenberg_matches_ed_and_python(version, monkeypatch):
    _require_julia(version)
    from dmrgpy.mpsjulialive import groundstate as jlgs
    rec = _Recorder(jlgs.Mainjl)
    monkeypatch.setattr(jlgs, "Mainjl", rec)
    n = 12
    energies = []
    for seed in (1, 2):
        _seed(seed)
        sc = _heisenberg(n, version, maxm=40, nsweeps=8, dmrg_method="CBE")
        energies.append(_quiet(sc.gs_energy))
    assert "get_gs_cbe" in rec.calls
    assert "get_gs_dmrg" not in rec.calls
    e_ed = sc.gs_energy(mode="ED")
    np.random.seed(1)
    sp = _heisenberg(n, "python", maxm=40, nsweeps=8, dmrg_method="CBE")
    e_py = _quiet(sp.gs_energy)
    for e in energies:
        assert e == pytest.approx(e_ed, abs=1e-8)
        assert e == pytest.approx(e_py, abs=1e-8)
    # the stored state is the solve's, at the reported energy
    assert sc.vev(sc.hamiltonian).real == pytest.approx(e_ed, abs=1e-8)


@pytest.mark.parametrize("version", JULIA)
def test_julia_cbe_spinless_fermions_match_ed(version):
    _require_julia(version)
    n = 8
    _seed(3)
    fc = fermionchain.Fermionic_Chain(n, itensor_version=version, maxm=40,
                                      nsweeps=8, dmrg_method="CBE")
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i]
        h = h + 0.7 * fc.N[i] * fc.N[i + 1]
    for i in range(n):
        h = h + (0.2 * i - 0.5) * fc.N[i]
    fc.set_hamiltonian(h)
    e = _quiet(fc.gs_energy)
    assert e == pytest.approx(fc.gs_energy(mode="ED"), abs=1e-8)
    assert fc.vev(fc.N[0]).real == pytest.approx(
        fc.vev(fc.N[0], mode="ED").real, abs=1e-6)


@pytest.mark.parametrize("version", JULIA)
def test_julia_cbe_grows_the_bonds_of_a_product_state(version):
    """From the Neel state every bond is 1, and one-site DMRG alone could
    only rotate each site. The expansion is what lets the bonds grow."""
    _require_julia(version)
    from dmrgpy.mpsjulialive.juliasession import Main
    from dmrgpy.mpsjulialive.mps import MPS
    n = 10
    sc = _heisenberg(n, version, maxm=32, nsweeps=6, dmrg_method="CBE")
    neel = Main.seval('sites -> MPS(sites, [isodd(i) ? "Up" : "Dn" '
                      'for i in 1:length(sites)])')(sc.jlsites)
    wf = MPS(neel, MBO=sc)
    assert max(_links(wf)) == 1
    e = _quiet(lambda: sc.gs_energy(wf0=wf))
    assert max(_links(sc.wf0)) > 1
    assert e == pytest.approx(sc.gs_energy(mode="ED"), abs=1e-8)


def _neel_states(sc, n):
    """The Neel state as a Julia product MPS (all bonds 1) and as an ED
    state, the ground state of a staggered field."""
    from dmrgpy.mpsjulialive.juliasession import Main
    from dmrgpy.mpsjulialive.mps import MPS
    h0 = 0
    for i in range(n):
        h0 = h0 + (-1) ** i * sc.Sz[i]
    sc.set_hamiltonian(h0)
    wfed = sc.get_gs(mode="ED")
    # site 1 (Sz[0]) down, the ground state of the field above
    neel = Main.seval('sites -> MPS(sites, [isodd(i) ? "Dn" : "Up" '
                      'for i in 1:length(sites)])')(sc.jlsites)
    return MPS(neel, MBO=sc), wfed


def _quench(method, n=10, nt=40, dt=0.1):
    sc = _heisenberg(n, "julia_live", maxm=32, nsweeps=6, tevol_method=method)
    h = sc.hamiltonian
    wf, wfed = _neel_states(sc, n)
    sc.set_hamiltonian(h + 0.2 * sc.Sz[1])
    op = sc.Sz[n // 2] + sc.Sz[n // 2 + 1] * 0.5
    _, z, wft = _quiet(lambda: timedependent.evolve_and_measure(
        sc, operator=op, wf=wf, nt=nt, dt=dt, return_wf=True))
    _, zed = timedependent.evolve_and_measure(sc, mode="ED", operator=op,
                                              wf=wfed, nt=nt, dt=dt)
    err = np.max(np.abs(np.real(np.asarray(z)) - np.real(np.asarray(zed))))
    return err, wft


@pytest.mark.parametrize("version", JULIA)
def test_julia_tdvp_cbe_quench_from_a_product_state_matches_ed(version, monkeypatch):
    """Bond dimension 1 at t=0, so one-site TDVP alone would never leave
    the product manifold; CBE-TDVP follows ED as two-site TDVP does."""
    _require_julia(version)
    from dmrgpy.mpsjulialive import juliasession
    rec = _Recorder(juliasession.Main)
    monkeypatch.setattr(juliasession, "Main", rec)
    err_cbe, wft = _quench("TDVP_CBE")
    assert "evolve_and_measure_tdvp_cbe" in rec.calls
    monkeypatch.undo()
    err_2s, _ = _quench("TDVP")
    print("TDVP_CBE error %.3e, TDVP error %.3e" % (err_cbe, err_2s))
    assert max(_links(wft)) > 1
    assert err_cbe < 5e-5
    assert err_cbe < 2 * err_2s + 1e-6


@pytest.mark.parametrize("version", JULIA)
def test_julia_tdvp_cbe_complex_time_step_matches_exact(version):
    """The one-step route TDZ takes (advance_complex_time_step), with a
    complex step, against exp(-i dz H) applied exactly."""
    _require_julia(version)
    import scipy.linalg as sla
    from dmrgpy.multioperator import MO2matrix
    from dmrgpy.mpsjulialive import timedependent as jtd
    from dmrgpy.mpsjulialive.mpo import MPO
    n = 6
    sc = _heisenberg(n, version, maxm=64, nsweeps=6, tevol_method="TDVP_CBE")
    h = sc.hamiltonian
    wf, wfed = _neel_states(sc, n)
    sc.set_hamiltonian(h)
    dz = 0.05 * np.exp(-0.3j)
    Hop = MPO(h, MBO=sc)
    for _ in range(4):
        wf = jtd.advance_complex_time_step(sc, Hop, wf, dz)
    Hm = np.array(MO2matrix(h, sc.get_ED_obj()).todense())
    v = sla.expm(-1j * 4 * dz * Hm) @ np.asarray(wfed.v).reshape(-1)
    for i in range(n):
        Sz = np.array(MO2matrix(sc.Sz[i], sc.get_ED_obj()).todense())
        exact = (v.conj() @ Sz @ v).real / (v.conj() @ v).real
        got = (wf.aMb(sc.Sz[i], wf) / wf.dot(wf)).real
        assert got == pytest.approx(exact, abs=1e-6)


@pytest.mark.parametrize("version", JULIA)
def test_julia_tdvp_cbe_correlator_matches_ed(version):
    """The quench correlator <GS|Sz e^{-i(H-E0)t} Sz|GS> (evolution_DC,
    the route of submode="TD"), against ED."""
    _require_julia(version)
    n = 8
    _seed(4)
    sc = _heisenberg(n, version, maxm=32, nsweeps=8, tevol_method="TDVP_CBE")
    _, c = _quiet(lambda: timedependent.evolution_DC(
        sc, name=(sc.Sz[3], sc.Sz[3]), nt=30, dt=0.1))
    _, ced = timedependent.evolution_DC(sc, mode="ED",
                                        name=(sc.Sz[3], sc.Sz[3]), nt=30, dt=0.1)
    err = np.max(np.abs(np.asarray(c) - np.asarray(ced)))
    print("TDVP_CBE correlator error %.3e" % err)
    assert err < 1e-4


def test_julia_cbe_checks_accept_julia_live():
    """Without starting Julia: both checks let julia_live through, and v2
    is still refused."""
    stub = types.SimpleNamespace(dmrg_method="CBE", cbe_delta=0.1,
                                 itensor_version="julia_live")
    groundstate.check_dmrg_method(stub)
    stub = types.SimpleNamespace(tevol_method="TDVP_CBE",
                                 itensor_version="julia_live")
    timedependent.check_tevol_method(stub)
    stub = types.SimpleNamespace(tevol_method="TDVP_CBE", itensor_version=2)
    with pytest.raises(NotImplementedError, match='"python"'):
        timedependent.check_tevol_method(stub)
