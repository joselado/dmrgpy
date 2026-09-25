"""Regression tests for the `session` cluster of the third 2026-09-24 hole
hunt (docs/audit_2026_09_24c_hole_hunt.md, findings 1 to 5, 8, 9 and 11).

   1. After set_gs() of a state off the ground manifold, DMRG KPM measured
      it from the solved E_0 while every other submode used the state's own
      energy (a rigid 0.3 shift on the chain below), and mode="ED" used
      three origins at once. Every submode on both modes now measures from
      the state's own energy; the KPM window stays on the band edges.
   2. mode="ED" submode="ED" read no state, neither set_gs() nor wf0=.
   3. set_gs() inside a sector followed directly by promote_to_dense()
      discarded the state; 4. after promote_to_dense() the next read
      re-swept any carried state that was not the sector ground state.
   5. After set_hamiltonian(H2, restart=False) every read but the public
      correlator answered for H1.
   8. submode="SECTOR" ignored set_gs().
   9. The first read after an injection ran an upper-edge solve and a
      discarded fluctuation; the lower edge no longer sweeps the state.
  11. A malformed KPM call paid a full ground-state solve before raising,
      and SECTOR paid one it never read.

Anchors are exact Lehmann sums and ED on chains of 3 and 4 sites.
"""

import io
import contextlib
import warnings

import numpy as np
import pytest

from dmrgpy import cppext, groundstate, spinchain
from dmrgpy.edtk import dynamics as eddynamics


def _backend(version):
    marks = ()
    if version in (2, 3):
        marks = pytest.mark.skipif(not cppext.available(version),
            reason="needs the compiled itensor_version=%d extension" % version)
    ident = "v%d" % version if version in (2, 3) else str(version)
    return pytest.param(version, id=ident, marks=marks)


BACKENDS = [_backend("python"), _backend(3)]
ALL = [_backend("python"), _backend(3), _backend(2)]


def _quiet(f):
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return f()


def _heis(sc, B=0.0):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(sc.ns):
        h = h + B*sc.Sz[i]
    return h


# ------------------------------------------------------ findings 1 and 2

N3, B3, DELTA = 3, 0.3, 0.05
ES = np.linspace(-0.8, 2.4, 641)


@pytest.fixture(scope="module")
def exact3():
    """The 3-site Heisenberg chain in Bz=0.3 (E = -1.15, -0.85, -0.15, 0.05,
    ...) and the exact (S+_0,S-_0) density of its first excited eigenstate
    |1> measured from E_1, with lines at -0.3, 0.7 and 1.2."""
    ed = spinchain.Spin_Chain(["S=1/2"]*N3)
    ed.set_hamiltonian(_heis(ed, B3))
    edo = ed.get_ED_obj()
    dense = lambda m: np.asarray(edo.MO2matrix(m).todense())
    E, V = np.linalg.eigh(dense(_heis(ed, B3)))
    Am, Bm = dense(ed.Sx[0] + 1j*ed.Sy[0]), dense(ed.Sx[0] - 1j*ed.Sy[0])
    st = V[:, 1]
    M = (st.conj() @ Am @ V)*(V.conj().T @ Bm @ st)
    own = sum(M[k]*DELTA/np.pi/((ES - (E[k] - E[1]))**2 + DELTA**2)
              for k in range(len(E)))
    return E, own


def _lines(y, thr=0.05):
    y = np.real(y); m = np.max(y)
    return [round(ES[k], 2) for k in range(1, len(y)-1)
            if y[k] >= y[k-1] and y[k] > y[k+1] and y[k] > thr*m]


def _pair(c):
    return (c.Sx[0] + 1j*c.Sy[0], c.Sx[0] - 1j*c.Sy[0])


@pytest.mark.parametrize("version", ALL)
@pytest.mark.parametrize("submode", ["KPM", "CVM", "ROOTN", "TD", "EX"])
def test_every_dmrg_submode_measures_a_set_state_from_its_own_energy(exact3, version, submode):
    E, own = exact3
    np.random.seed(2)
    sc = spinchain.Spin_Chain(["S=1/2"]*N3, itensor_version=version)
    sc.maxm, sc.nsweeps = 20, 12
    sc.set_hamiltonian(_heis(sc, B3))
    sc.gs_energy()
    ee, ww = _quiet(lambda: sc.get_excited_states(n=2))
    sc.set_gs(ww[1])
    kw = dict(delta=DELTA, es=ES)
    if submode == "EX": kw["nex"] = 8
    if submode == "TD": kw["dt"] = 0.05
    x, y = _quiet(lambda: sc.get_dynamical_correlator(submode=submode, name=_pair(sc), **kw))
    assert sc.gs_energy() == pytest.approx(E[1], abs=1e-6)
    assert _lines(y) == [-0.3, 0.7, 1.2]
    if submode != "KPM": # KPM's Jackson line is not a Lorentzian
        assert np.max(np.abs(y - own)) < 2e-2*np.max(np.abs(own))


@pytest.mark.parametrize("submode", ["KPM", "CVM", "INV", "ROOTN", "TD", "EX", "ED"])
def test_every_ed_submode_measures_a_set_state_from_its_own_energy(exact3, submode):
    """KPM, CVM, INV, ROOTN and TD measured it from the solved E_0, and
    submode="ED" read the ground state instead (2.834 off on a 2.835
    peak)."""
    E, own = exact3
    ed = spinchain.Spin_Chain(["S=1/2"]*N3)
    ed.set_hamiltonian(_heis(ed, B3))
    ed.get_gs(mode="ED")
    ee, ww = ed.get_excited_states(n=2, mode="ED")
    ed.set_gs(ww[1])
    kw = dict(delta=DELTA, es=ES)
    if submode == "EX": kw["nex"] = 8
    x, y = _quiet(lambda: ed.get_dynamical_correlator(mode="ED", submode=submode,
                                                    name=_pair(ed), **kw))
    assert _lines(y) == [-0.3, 0.7, 1.2]
    if submode in ("ED", "INV", "CVM", "EX"):
        assert np.max(np.abs(y - own)) < 1e-6*np.max(np.abs(own))


def test_ed_submode_ed_takes_an_explicit_state():
    """wf0= was computed two lines above the ED branch and never forwarded."""
    ed = spinchain.Spin_Chain(["S=1/2"]*N3)
    ed.set_hamiltonian(_heis(ed, B3))
    ee, ww = ed.get_excited_states(n=2, mode="ED")
    x, y = ed.get_dynamical_correlator(mode="ED", submode="ED", name=_pair(ed),
                                       wf0=ww[1], delta=DELTA, es=ES)
    assert _lines(y) == [-0.3, 0.7, 1.2]


def test_ed_submode_ed_default_is_the_manifold_average():
    """With no state set or passed, the dex manifold average stands,
    bit for bit, on a degenerate ground doublet (B=0)."""
    ed = spinchain.Spin_Chain(["S=1/2"]*N3)
    ed.set_hamiltonian(_heis(ed, 0.0))
    x, y = ed.get_dynamical_correlator(mode="ED", submode="ED", name=_pair(ed),
                                       delta=DELTA, es=ES)
    edo = ed.get_ED_obj()
    emu, vs = edo.get_diagonalized_hamiltonian()
    A = edo.MO2matrix(_pair(ed)[0]); B = edo.MO2matrix(_pair(ed)[1])
    _, ref = eddynamics.dynamical_correlator_ED(edo.get_hamiltonian(), A, B,
                                                 emu=emu, vs=vs, delta=DELTA, es=ES)
    assert np.array_equal(y, ref)


@pytest.mark.parametrize("version", BACKENDS)
def test_kpm_axis_at_a_solved_ground_state_is_unchanged(version):
    """The origin moved from the band edge emin to the state's e0, which
    are the same number at a solved ground state: the public curve equals
    the reconstruction on emin exactly."""
    from dmrgpy import kpmdmrg
    np.random.seed(3)
    sc = spinchain.Spin_Chain(["S=1/2"]*6, itensor_version=version)
    sc.maxm, sc.nsweeps = 20, 10
    sc.set_hamiltonian(_heis(sc, 0.2))
    x, y = sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]), delta=0.1, es=ES)
    mus, emin, emax, scale, n, d = sc.get_dynamical_correlator_moments(
        name=(sc.Sz[0], sc.Sz[0]), delta=0.1)
    _, yref = kpmdmrg.dynamical_correlator_from_moments(mus, emin, emax, scale, n, ES,
                                                        delta=d)
    assert emin == pytest.approx(sc.gs_energy(), abs=1e-12)
    assert np.max(np.abs(y - yref)) < 1e-10*np.max(np.abs(yref))


# ------------------------------------------------------ findings 3 and 4

def _sector_chain(version, seed):
    np.random.seed(seed)
    sc = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version=version)
    sc.maxm, sc.nsweeps = 20, 12
    sc.set_hamiltonian(_heis(sc, -0.9))
    sc.set_conserved_sector(Sz=0)
    sc.gs_energy()
    ee, ww = _quiet(lambda: sc.get_excited_states(n=2))
    x = (ww[1] + 0.6*ww[0]).normalize()
    ex = (sc.aMb(x, sc.hamiltonian, x)/sc.overlap(x, x)).real
    return sc, x, ex


def _kept(sc, x):
    xd = sc.promote_mps(x); wf = sc.wf0
    return abs(sc.overlap(xd, wf))**2/(abs(sc.overlap(xd, xd))*abs(sc.overlap(wf, wf)))


@pytest.mark.parametrize("version", BACKENDS)
@pytest.mark.parametrize("setter", ["set_gs", "set_initial_wf"])
@pytest.mark.parametrize("read_first", [False, True])
def test_promote_to_dense_keeps_a_set_state(version, setter, read_first):
    """x = normalize(|1> + 0.6|0>) in Sz=0, in a field that puts the global
    ground state at Sz_tot=1: finding 3 lost it with no read in between
    (overlap 0.0000), finding 4 with one (the global ground state -1.857107
    on "python", the sector's -1.616025 on v3)."""
    sc, x, ex = _sector_chain(version, 3)
    getattr(sc, setter)(x)
    if read_first: sc.gs_energy()
    sc.promote_to_dense()
    assert sc.gs_energy() == pytest.approx(ex, abs=1e-8)
    assert _kept(sc, x) == pytest.approx(1.0, abs=1e-8)
    assert sum(sc.vev(sc.Sz[i]).real for i in range(4)) == pytest.approx(0.0, abs=1e-8)


@pytest.mark.parametrize("version", BACKENDS)
def test_promote_to_dense_keeps_the_solved_sector_state(version):
    np.random.seed(4)
    sc = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version=version)
    sc.maxm, sc.nsweeps = 20, 12
    sc.set_hamiltonian(_heis(sc, -0.9))
    sc.set_conserved_sector(Sz=0)
    e = sc.gs_energy()
    sc.promote_to_dense()
    assert sc.gs_energy() == pytest.approx(e, abs=1e-8)
    assert e == pytest.approx(-1.616025, abs=1e-6)


# ------------------------------------------------------------ finding 5

@pytest.mark.parametrize("version", ALL)
def test_restart_false_answers_for_the_new_hamiltonian_everywhere(version):
    """4-site Heisenberg plus 0.8*Sx_0, which does not commute with H1:
    gs_energy() read -1.616025 against -1.780099 and <Sx_0> 0 against
    -0.349221, and get_excited() gave H2's Ritz values in H1's states."""
    np.random.seed(5)
    sc = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version=version)
    sc.maxm, sc.nsweeps = 20, 12
    h1 = _heis(sc)
    sc.set_hamiltonian(h1)
    sc.gs_energy()
    sc.set_hamiltonian(h1 + 0.8*sc.Sx[0], restart=False)
    ev = np.linalg.eigvalsh(sc.get_ED_obj().get_hamiltonian().toarray())
    assert sc.gs_energy() == pytest.approx(ev[0], abs=1e-6)
    assert sc.gs_energy() == pytest.approx(-1.780099, abs=1e-6)
    assert sc.vev(sc.Sx[0]).real == pytest.approx(-0.349221, abs=1e-5)
    es = _quiet(lambda: sc.get_excited(n=2))
    assert np.real(es) == pytest.approx(ev[:2], abs=1e-5)


@pytest.mark.parametrize("version", BACKENDS)
def test_restart_false_direct_kpm_moments_are_the_new_hamiltonians(version):
    np.random.seed(6)
    sc = spinchain.Spin_Chain(["S=1/2"]*6, itensor_version=version)
    sc.maxm, sc.nsweeps = 20, 12
    h1 = _heis(sc)
    sc.set_hamiltonian(h1)
    sc.gs_energy()
    h2 = h1 + sum((1.0 if i % 2 else -1.0)*sc.Sz[i] for i in range(6))
    sc.set_hamiltonian(h2, restart=False)
    mus, emin, emax, *_ = sc.get_dynamical_correlator_moments(
        name=(sc.Sz[0], sc.Sz[0]), delta=0.2)
    fresh = spinchain.Spin_Chain(["S=1/2"]*6, itensor_version=version)
    fresh.maxm, fresh.nsweeps = 20, 12
    fresh.set_hamiltonian(h2)
    mus2, emin2, *_ = fresh.get_dynamical_correlator_moments(
        name=(fresh.Sz[0], fresh.Sz[0]), delta=0.2)
    assert emin == pytest.approx(emin2, abs=1e-8)
    assert len(mus) == len(mus2)
    assert np.max(np.abs(np.array(mus) - np.array(mus2))) < 1e-6


# ------------------------------------------------------------ finding 8

@pytest.mark.parametrize("version", BACKENDS)
def test_sector_refuses_a_set_state_it_does_not_measure(version):
    """After set_gs(|2>), the upper doublet's Sz=-1/2 member, SECTOR
    returned the global ground state's spectrum, 1.0485 off on a 1.0616
    peak."""
    np.random.seed(7)
    sc = spinchain.Spin_Chain(["S=1/2"]*N3, itensor_version=version)
    sc.maxm, sc.nsweeps = 20, 12
    sc.set_hamiltonian(_heis(sc, B3))
    sc.gs_energy()
    ee, ww = _quiet(lambda: sc.get_excited_states(n=4))
    k = int(np.argmin(np.abs(np.real(ee) - (-0.15))))
    sc.set_gs(ww[k])
    with pytest.raises(NotImplementedError, match="not that state"):
        _quiet(lambda: sc.get_dynamical_correlator(submode="SECTOR",
               name=(sc.Sz[0], sc.Sz[0]), delta=DELTA, es=ES, nex=4))


@pytest.mark.parametrize("version", BACKENDS)
def test_sector_measures_a_set_state_that_is_its_sectors_ground_state(exact3, version):
    """|1>, the Sz=+1/2 member, is the lowest state of its own sector:
    SECTOR now reads its charge from the chain's state, not from the
    clone's unconstrained re-solve, and returns its curve."""
    E, own = exact3
    np.random.seed(8)
    sc = spinchain.Spin_Chain(["S=1/2"]*N3, itensor_version=version)
    sc.maxm, sc.nsweeps = 20, 12
    sc.set_hamiltonian(_heis(sc, B3))
    sc.gs_energy()
    ee, ww = _quiet(lambda: sc.get_excited_states(n=2))
    sc.set_gs(ww[1])
    x, y = _quiet(lambda: sc.get_dynamical_correlator(submode="SECTOR",
                  name=_pair(sc), delta=DELTA, es=ES, nex=4))
    assert np.max(np.abs(y - own)) < 1e-6*np.max(np.abs(own))


# ------------------------------------------------------------ finding 9

@pytest.mark.parametrize("version", BACKENDS)
def test_an_injected_read_fills_no_band_edge(version, monkeypatch):
    """The first read after set_gs() ran session.excited_states(1): an
    upper-edge solve and a discarded fluctuation, 7.1 s on 24 sites of
    "python" for a 0.04 s <x|H|x>."""
    np.random.seed(9)
    sc = spinchain.Spin_Chain(["S=1/2"]*6, itensor_version=version)
    sc.maxm, sc.nsweeps = 20, 8
    sc.set_hamiltonian(_heis(sc))
    x = sc.random_mps()
    calls = []
    sess = sc._session
    orig = sess.excited_states
    monkeypatch.setattr(sess, "excited_states", lambda *a, **k: calls.append(a) or orig(*a, **k),
                        raising=False) if version == "python" else None
    sc.set_gs(x)
    e = sc.gs_energy()
    assert e == pytest.approx((sc.aMb(x, sc.hamiltonian, x)/sc.overlap(x, x)).real, abs=1e-12)
    if version == "python":
        assert calls == []
        assert sess._bandwidth_max is None


@pytest.mark.parametrize("version", BACKENDS)
def test_the_first_kpm_call_after_an_injection_does_not_sweep_it(version):
    """What the removed pre-fill was for: the lower band edge is now filled
    from a fresh start, so the injected state survives a KPM call."""
    np.random.seed(10)
    sc = spinchain.Spin_Chain(["S=1/2"]*3, itensor_version=version)
    sc.maxm, sc.nsweeps = 20, 10
    sc.set_hamiltonian(_heis(sc))
    x = sc.random_mps()
    sc.set_gs(x)
    sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]), delta=0.2, es=ES)
    wf = sc.get_gs()
    assert abs(x.dot(wf))**2/(abs(x.dot(x))*abs(wf.dot(wf))) == pytest.approx(1.0, abs=1e-10)
    # and the window's bottom is the Hamiltonian's ground energy
    mus, emin, *_ = sc.get_dynamical_correlator_moments(name=(sc.Sz[0], sc.Sz[0]), delta=0.2)
    ev = np.linalg.eigvalsh(sc.get_ED_obj().get_hamiltonian().toarray())
    assert emin == pytest.approx(ev[0], abs=1e-8)


# ------------------------------------------------------------ finding 11

def _count_solves(monkeypatch):
    calls = []
    orig = groundstate.gs_energy_single

    def spy(*a, **k):
        calls.append(1)
        return orig(*a, **k)
    monkeypatch.setattr(groundstate, "gs_energy_single", spy)
    return calls


@pytest.mark.parametrize("version", ALL)
@pytest.mark.parametrize("bad", ["n_scale", "keyword", "delta"])
def test_malformed_kpm_call_raises_before_the_ground_state(version, bad, monkeypatch):
    sc = spinchain.Spin_Chain(["S=1/2"]*6, itensor_version=version)
    sc.set_hamiltonian(_heis(sc))
    calls = _count_solves(monkeypatch)
    kw = dict(name=(sc.Sz[0], sc.Sz[0]))
    if bad == "n_scale": sc.kpm_n_scale = 1.5; err = TypeError
    elif bad == "keyword": kw["deltaa"] = 0.1; err = TypeError
    else: kw["delta"] = -0.1; err = ValueError
    with pytest.raises(err):
        sc.get_dynamical_correlator(**kw)
    assert calls == [] and not sc.computed_gs


@pytest.mark.skipif(not cppext.available(2), reason="needs the compiled v2 extension")
def test_v2_energy_truncation_raises_before_the_ground_state(monkeypatch):
    sc = spinchain.Spin_Chain(["S=1/2"]*6, itensor_version=2)
    sc.set_hamiltonian(_heis(sc))
    sc.kpm_energy_truncate = True
    calls = _count_solves(monkeypatch)
    with pytest.raises(NotImplementedError):
        sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]))
    assert calls == [] and not sc.computed_gs


@pytest.mark.parametrize("version", BACKENDS)
def test_sector_makes_no_solve_on_the_callers_chain(version):
    np.random.seed(11)
    sc = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version=version)
    sc.maxm, sc.nsweeps = 20, 10
    sc.set_hamiltonian(_heis(sc))
    _quiet(lambda: sc.get_dynamical_correlator(submode="SECTOR",
           name=(sc.Sx[0] - 1j*sc.Sy[0], sc.Sx[0] + 1j*sc.Sy[0]), delta=0.2, es=ES, nex=4))
    assert not sc.computed_gs


@pytest.mark.parametrize("version", BACKENDS)
def test_kpm_under_energy_truncation_measures_a_set_state_from_its_own_energy(version):
    """The one route where both numbers finding 1 separated enter: the
    gs-anchored window is built on the Hamiltonian's E_0 (a fresh-start
    minimum_energy since finding 9) and the axis on the set state's e0.
    The lines sit at -0.3 and 0.7, not at 0.0 and 1.0; the 1.2 line is
    absent under the truncation on both backends, which was not examined
    (docs/known_issue_kpm_energy_truncation_window.md's territory)."""
    np.random.seed(2)
    sc = spinchain.Spin_Chain(["S=1/2"]*N3, itensor_version=version)
    sc.maxm, sc.nsweeps = 20, 12
    sc.set_hamiltonian(_heis(sc, B3))
    sc.gs_energy()
    ee, ww = _quiet(lambda: sc.get_excited_states(n=2))
    sc.kpm_energy_truncate = True
    sc.set_gs(ww[1])
    x, y = _quiet(lambda: sc.get_dynamical_correlator(name=_pair(sc), delta=DELTA, es=ES))
    got = _lines(y)
    assert -0.3 in got and 0.7 in got
    assert 0.0 not in got and 1.0 not in got
