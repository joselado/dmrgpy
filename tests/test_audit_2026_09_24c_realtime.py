"""Regression tests for the `realtime` cluster of the third 2026-09-24 hole
hunt (docs/audit_2026_09_24c_hole_hunt.md, findings 14 to 17).

  14. On itensor_version=3, tevol_method="TDVP_GSE" lost its Krylov
      expansion at the left edge before the first one-site update when
      site 0 was pinned to one local basis state, which every ladder
      operator or projector on site 0 produces, so site 0 stayed frozen for
      the whole trajectory.
  15. evolve_and_measure(mode="DMRG") and evolution_ABA(mode="DMRG") took a
      **kwargs nothing read, so DT=0.2 ran silently at dt=1e-2.
  16. On mode="ED" the same two raised TypeError on h=, which DMRG honours,
      and defaulted to nt=100 against DMRG's 1000.
  17. The ED propagator was RK45 at scipy's default tolerances, with an
      error set by dt times the absolute energy: a constant +20 added to H
      moved evolve_and_measure(mode="ED") from 2.1e-7 to 5.6e-4 off exact
      at dt=0.1. It is expm_multiply now.

Anchors: closed forms (Larmor precession) and a dense e^{-iHt} from the ED
matrix.
"""

import numpy as np
import pytest
import scipy.linalg as dlg

from dmrgpy import cppext, fermionchain, multioperator, spinchain, timedependent


def _backend(version):
    marks = ()
    if version in (2, 3):
        marks = pytest.mark.skipif(not cppext.available(version),
            reason="needs the compiled itensor_version=%d extension" % version)
    ident = "v%d" % version if version in (2, 3) else str(version)
    return pytest.param(version, id=ident, marks=marks)


BACKENDS = [_backend("python"), _backend(3), _backend(2)]


def _larmor(version):
    """3 sites, ground state of -sum Sx (the +x state), to be evolved under
    H2 = sum Sz, so <Sx_0>(t) = cos(t)/2 exactly."""
    sc = spinchain.Spin_Chain(["S=1/2"]*3, itensor_version=version)
    h1 = -sum(sc.Sx[i] for i in range(3))
    h2 = sum(sc.Sz[i] for i in range(3))
    sc.set_hamiltonian(h1)
    sc.maxm, sc.nsweeps = 8, 10
    np.random.seed(0)
    sc.gs_energy()
    return sc, h1, h2


# ------------------------------------------------------------ finding 15

@pytest.mark.parametrize("version", BACKENDS)
@pytest.mark.parametrize("entry", ["evolve_and_measure", "evolution_ABA"])
def test_misspelled_keyword_raises_on_dmrg(version, entry):
    sc, h1, h2 = _larmor(version)
    call = lambda **k: getattr(timedependent, entry)(sc, **k)
    kw = dict(operator=sc.Sx[0]) if entry == "evolve_and_measure" else dict(B=sc.Sx[0])
    with pytest.raises(TypeError, match="DT"):
        call(mode="DMRG", nt=20, DT=0.2, **kw)
    ts, cs = call(mode="DMRG", nt=5, dt=0.2, **kw)
    assert ts[-1] == pytest.approx(0.8)


# ------------------------------------------------------------ finding 16

@pytest.mark.parametrize("entry", ["evolve_and_measure", "evolution_ABA"])
def test_h_keyword_on_ed_evolves_under_it(entry):
    """h=H2 used to raise 'got multiple values for argument h' on ED."""
    sc, h1, h2 = _larmor("python")
    call = lambda **k: getattr(timedependent, entry)(sc, **k)
    kw = dict(operator=sc.Sx[0]) if entry == "evolve_and_measure" else dict(B=sc.Sx[0])
    ts, ed = call(mode="ED", h=h2, nt=5, dt=0.2, **kw)
    _, dm = call(mode="DMRG", h=h2, nt=5, dt=0.2, **kw)
    assert np.max(np.abs(ed - np.cos(ts)/2)) < 1e-12
    assert np.max(np.abs(dm - ed)) < 1e-10


@pytest.mark.parametrize("entry", ["evolve_and_measure", "evolution_ABA"])
def test_default_time_grid_is_the_same_on_both_modes(entry):
    sc, h1, h2 = _larmor("python")
    call = lambda **k: getattr(timedependent, entry)(sc, **k)
    kw = dict(operator=sc.Sx[0]) if entry == "evolve_and_measure" else dict(B=sc.Sx[0])
    ts_ed, _ = call(mode="ED", **kw)
    ts_dm, _ = call(mode="DMRG", **kw)
    assert len(ts_ed) == len(ts_dm) == 1000
    assert np.array_equal(ts_ed, ts_dm)


# ------------------------------------------------------------ finding 17

def _xxz_chain(n=8):
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version="python")
    h = sum(3.0*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]) + 3.0*sc.Sz[i]*sc.Sz[i+1]
            for i in range(n-1)) + 3.0*sum(sc.Sz[i] for i in range(0, n, 2))
    sc.set_hamiltonian(h)
    return sc, h


def test_ed_evolution_does_not_depend_on_the_energy_origin():
    """A constant added to H is a global phase; RK45 at its default
    tolerances resolved it as a frequency, 5.6e-4 off at +20, dt=0.1."""
    sc, h = _xxz_chain(8)
    _, a = timedependent.evolve_and_measure(sc, mode="ED", operator=sc.Sz[2], nt=40, dt=0.1)
    _, b = timedependent.evolve_and_measure(sc, mode="ED", operator=sc.Sz[2], nt=40, dt=0.1,
                                 h=h + 20.0*multioperator.identity())
    assert np.max(np.abs(a - b)) < 1e-12


def test_ed_evolution_matches_the_dense_propagator():
    """evolution_ABA(mode="ED") against a hand-written e^{-iHt}: 8.46e-4 off
    on a 0.068 scale at dt=0.2 before."""
    sc, h = _xxz_chain(8)
    ed = sc.get_ED_obj()
    H = ed.MO2matrix(h).toarray()
    A = ed.MO2matrix(sc.Sx[0]).toarray()
    B = ed.MO2matrix(sc.Sz[3]).toarray()
    w, v = np.linalg.eigh(H)
    gs = v[:, 0]
    nt, dt = 20, 0.2
    ts, cs = timedependent.evolution_ABA(sc, mode="ED", A=sc.Sx[0], B=sc.Sz[3], nt=nt, dt=dt)
    psi0 = A @ gs
    ref = []
    for t in ts:
        psi = v @ (np.exp(-1j*w*t)*(v.conj().T @ psi0))
        ref.append(psi.conj() @ B @ psi)
    assert np.max(np.abs(cs - np.array(ref))) < 1e-12


# ------------------------------------------------------------ finding 14

N14, NT14, DT14 = 6, 40, 0.05
TS14 = DT14*np.arange(NT14)
V3 = pytest.mark.skipif(not cppext.available(3), reason="needs the compiled v3 extension")


def _ferm14(version):
    fc = fermionchain.Fermionic_Chain(N14, itensor_version=version)
    fc.maxm, fc.nsweeps = 40, 20
    fc.tevol_method = "TDVP_GSE"
    h = 0
    for j in range(N14 - 1):
        h = h - fc.Cdag[j+1]*fc.C[j] - fc.Cdag[j]*fc.C[j+1]
    h = h + 0.3*fc.N[0] + 0.5*fc.N[2]*fc.N[3] - 0.8*sum(fc.N[j] for j in range(N14))
    fc.set_hamiltonian(h)
    return fc


def _spin14(version):
    sc = spinchain.Spin_Chain([2]*N14, itensor_version=version)
    sc.maxm, sc.nsweeps = 40, 20
    sc.tevol_method = "TDVP_GSE"
    h = sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
            for i in range(N14 - 1)) + 0.2*sc.Sz[0]
    sc.set_hamiltonian(h)
    return sc


def _exact_aba(ch, A, B):
    ed = ch.get_ED_obj()
    H = ed.get_operator(ch.hamiltonian).toarray()
    w, v = np.linalg.eigh(H)
    psi0 = ed.get_operator(A).toarray() @ v[:, 0]
    Bm = ed.get_operator(B).toarray()
    c = np.conj(v).T @ psi0
    return np.array([np.vdot(v @ (np.exp(-1j*w*t)*c), Bm @ (v @ (np.exp(-1j*w*t)*c)))
                     for t in TS14])


STARTS14 = [
    ("fermion", "C_0", lambda c: c.C[0], lambda c: c.N[0]),
    ("fermion", "1-N_0", lambda c: 1 - c.N[0], lambda c: c.N[0]),
    ("fermion", "N_0", lambda c: c.N[0], lambda c: c.N[0]),
    ("spin", "Sp_0", lambda c: c.Sx[0] + 1j*c.Sy[0], lambda c: c.Sz[0]),
    ("spin", "Sm_0", lambda c: c.Sx[0] - 1j*c.Sy[0], lambda c: c.Sz[0]),
    ("spin", "1/2+Sz_0", lambda c: 0.5 + c.Sz[0], lambda c: c.Sz[0]),
]


@V3
@pytest.mark.parametrize("kind,label,mkA,mkB", STARTS14, ids=[s[1] for s in STARTS14])
def test_v3_tdvp_gse_moves_a_site_pinned_at_the_left_edge(kind, label, mkA, mkB):
    """Site 0 of the evolved state stayed frozen on v3 (C_0: 4.74e-01 off,
    site-0 entropy 0 against 0.6330), because the Cutoff=0 sweeps after
    addBasis discarded the exactly-zero-weight directions GSE had added at
    the edge. "python" TDVP_GSE, the reference, was exact throughout."""
    ch = (_ferm14 if kind == "fermion" else _spin14)(3)
    A, B = mkA(ch), mkB(ch)
    ref = _exact_aba(ch, A, B)
    _t, y = timedependent.evolution_ABA(ch, A=A, B=B, mode="DMRG", nt=NT14, dt=DT14)
    assert np.ptp(np.real(ref)) > 1e-2 # site 0 does move
    assert np.max(np.abs(np.asarray(y) - ref)) < 1e-5


@V3
def test_v3_tdvp_gse_td_spectrum_at_the_default_sites():
    """The public (Cdag_0,C_0) TD spectrum on v3 TDVP_GSE was 1.735 times its
    own peak off (a height of 0.0808 against 0.4756); it now agrees with
    v3's own two-site TDVP."""
    es = np.linspace(-3.0, 3.0, 301)
    out = {}
    for method in ("TDVP_GSE", "TDVP"):
        fc = _ferm14(3)
        fc.tevol_method = method
        np.random.seed(0)
        fc.gs_energy()
        _, out[method] = fc.get_dynamical_correlator(submode="TD", predict=False,
            name=(fc.Cdag[0], fc.C[0]), delta=0.2, es=es, dt=0.1)
    peak = np.max(np.abs(out["TDVP"]))
    assert np.max(np.abs(out["TDVP_GSE"] - out["TDVP"])) < 1e-3*peak
