"""Regression tests for the `session` cluster of the 2026-09-25 fixes of the
open items the 2026-09-24 records left behind, all at the edges of the
injected-state machinery (groundstate.mark_injected and its readers).

  julia-warm-start. On julia_live, set_initial_wf() and
      set_initial_wf_guess() never reached the solver, which started from
      its own random state, and set_gs() raised AttributeError on the Julia
      MPS. julia_live now follows gs_energy_single()'s contract: a set state
      is taken unswept with e0 = <x|H|x>, a guess is swept from, the
      caller's state is not moved, and KPM measures a set state from its
      own energy on the band-edge window.
  generalized-cache. gs_energy_generalized() re-sent H behind the send
      cache, so on a chain whose cache was empty the next correlator
      re-solved a plain ground state over the generalized one; on a
      non-Hermitian H neither NH solve sent H at all, so the same happened
      whether or not the chain was solved first, and a plain NH solve was
      solved again by the first correlator after it.
  thermal-bypass. Thermal_Spin_Chain.get_gs() assigned MBChain.wf0 and
      MBChain.hamiltonian directly, so a correlator re-solved over the
      annealed state and vev() on mode="ED" read the singlet state.
  nh-injected. After set_gs(x) the non-Hermitian KPM paired x with the
      left state of an earlier NH-DMRG solve, or raised AttributeError, and
      gs_energy(wf0=x) dropped x and re-solved.
  maxde-per-site. maxde is a fluctuation per site while
      gs_energy_fluctuation() is the total, which no docstring said.

Anchors are exact eigenstates, <x|O|x> of the state set, ED, and exact
Boltzmann averages on chains of 3 to 10 sites.
"""

import io
import contextlib
import re
import warnings

import numpy as np
import pytest

from dmrgpy import cppext, groundstate, spinchain, thermal
from dmrgpy.manybodychain import Many_Body_Chain

from _helpers import julia_available


def _backend(version):
    marks = ()
    if version in (2, 3):
        marks = pytest.mark.skipif(not cppext.available(version),
            reason="needs the compiled itensor_version=%d extension" % version)
    ident = "v%d" % version if version in (2, 3) else str(version)
    return pytest.param(version, id=ident, marks=marks)


BACKENDS = [_backend("python"), _backend(3)]


def _require_julia(version):
    """The julia_live check, deferred to the test body so that collecting
    this file (and -k "not julia_live") never boots a Julia session"""
    if version == "julia_live":
        ok, reason = julia_available()
        if not ok:
            pytest.skip("requires a working juliacall/Julia toolchain: %s"
                        % reason)


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


def _chain(version, n=3, B=0.0, maxm=10, nsweeps=10):
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=version)
    sc.set_hamiltonian(_heis(sc, B))
    sc.maxm, sc.nsweeps = maxm, nsweeps
    return sc


def _energy(sc, wf):
    return float(np.real(wf.aMb(sc.hamiltonian, wf)/wf.dot(wf)))


def _expect(wf, op):
    return complex(wf.dot(op*wf)/wf.dot(wf))


def _fidelity(a, b):
    return abs(a.dot(b))**2/abs(a.dot(a)*b.dot(b))


def _normalized(wf):
    return wf*(1/np.sqrt(wf.dot(wf).real))


def _nh_chain(version):
    sc = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version=version)
    sc.set_hamiltonian(_heis(sc) + 0.3j*sc.Sz[0] + 0.2*sc.Sx[1])
    sc.maxm, sc.nsweeps = 20, 10
    return sc


NHKPM = dict(submode="KPM", es=np.linspace(0.0, 4.0, 5), delta=0.3,
             E_max=10.0, n=50)


# ------------------------------------------------------- julia-warm-start

@pytest.mark.parametrize("version", BACKENDS + [_backend("julia_live")])
def test_warm_start_setters_reach_the_solver(version):
    """The two exact members of the 3-site ground doublet, built from the
    solved state as (Sz_tot -+ 1/2)|s>: a guess is swept from and stays on
    its member, a set state is taken as it is, and a state that is not an
    eigenstate keeps its own energy and expectation values. On julia_live
    the guess used to land at |<target|gs>|^2 = 0.0008 to 0.68, and set_gs
    raised AttributeError."""
    _require_julia(version)
    np.random.seed(1)
    sc = _chain(version)
    _quiet(sc.gs_energy)
    s = sc.get_gs().copy()
    Szt = sc.Sz[0] + sc.Sz[1] + sc.Sz[2]
    up = _normalized((Szt + 0.5)*s)
    dn = _normalized((0.5 - 1*Szt)*s)
    target, other = (dn, up) if np.real(_expect(s, Szt)) > 0 else (up, dn)
    sc.set_initial_wf_guess(target)
    _quiet(sc.gs_energy)
    assert _fidelity(sc.get_gs(), target) == pytest.approx(1.0, abs=1e-8)
    sc.set_initial_wf(other)
    assert _quiet(sc.gs_energy) == pytest.approx(_energy(sc, other), abs=1e-10)
    assert _fidelity(sc.get_gs(), other) == pytest.approx(1.0, abs=1e-10)
    x = _normalized(s + 0.4*(sc.Sx[0]*s)) # not an eigenstate
    e_x = _energy(sc, x)
    sc.set_gs(x)
    assert _quiet(sc.gs_energy) == pytest.approx(e_x, abs=1e-10)
    assert e_x > -1.0 + 1e-3 # so the solved -1 cannot pass for it
    zz = sc.Sz[0]*sc.Sz[1]
    assert np.real(sc.vev(zz)) == pytest.approx(np.real(_expect(x, zz)), abs=1e-10)
    # a sweep from the caller's state does not move it
    keep = x.copy()
    sc.set_initial_wf_guess(x)
    assert _quiet(sc.gs_energy) < e_x - 1e-3
    assert _energy(sc, x) == pytest.approx(e_x, abs=1e-12)
    assert _fidelity(x, keep) == pytest.approx(1.0, abs=1e-12)
    # the explicit warm start on a solved chain (julia_live raised TypeError)
    _quiet(lambda: sc.gs_energy(wf0=target))
    assert _fidelity(sc.get_gs(), target) == pytest.approx(1.0, abs=1e-8)


@pytest.mark.parametrize("version", [_backend("julia_live")])
def test_julia_kpm_measures_a_set_state_from_its_own_energy(version):
    """The 3-site chain in Bz=0.3 (E = -1.15, -0.85, ...) after set_gs of its
    first excited eigenstate |1>, the Sz=+1/2 member S+_tot|0>: the
    (S+_0,S-_0) lines sit at E_n - E_1 = -0.3, 0.7 and 1.2, the own-origin
    lines the session backends return (2026-09-24c audit, finding 1), and
    the curve is mode="ED"'s, which needs the window on the band edges:
    anchored on E_1 instead it is 1.28 of the peak off, against 8.2e-4."""
    _require_julia(version)
    ES = np.linspace(-0.8, 2.4, 641)
    KW = dict(submode="KPM", delta=0.05, es=ES)
    ed = spinchain.Spin_Chain(["S=1/2"]*3)
    ed.set_hamiltonian(_heis(ed, 0.3))
    ed.get_gs(mode="ED")
    ee, ww = ed.get_excited_states(n=2, mode="ED")
    ed.set_gs(ww[1])
    _, yed = _quiet(lambda: ed.get_dynamical_correlator(mode="ED",
            name=(ed.Sx[0] + 1j*ed.Sy[0], ed.Sx[0] - 1j*ed.Sy[0]), **KW))
    sc = _chain(version, B=0.3, maxm=20, nsweeps=12)
    _quiet(sc.gs_energy)
    Sp = sc.Sx[0] + 1j*sc.Sy[0]
    Spt = Sp + (sc.Sx[1] + 1j*sc.Sy[1]) + (sc.Sx[2] + 1j*sc.Sy[2])
    one = _normalized(Spt*sc.get_gs())
    assert _energy(sc, one) == pytest.approx(-0.85, abs=1e-6)
    sc.set_gs(one)
    x, y = _quiet(lambda: sc.get_dynamical_correlator(
        name=(Sp, sc.Sx[0] - 1j*sc.Sy[0]), **KW))
    yr = np.real(y); m = np.max(yr)
    lines = [round(ES[k], 2) for k in range(1, len(yr)-1)
             if yr[k] >= yr[k-1] and yr[k] > yr[k+1] and yr[k] > 0.05*m]
    assert lines == [-0.3, 0.7, 1.2]
    assert np.max(np.abs(y - yed)) < 1e-2*np.max(np.abs(yed))
    assert sc.gs_energy() == pytest.approx(-0.85, abs=1e-6)
    assert _fidelity(sc.get_gs(), one) == pytest.approx(1.0, abs=1e-10)


# ------------------------------------------------------ generalized-cache

@pytest.mark.parametrize("version", BACKENDS)
@pytest.mark.parametrize("solved_first", [False, True], ids=["fresh", "solved"])
def test_correlator_after_gs_energy_generalized_reads_the_generalized_state(
        version, solved_first):
    """On a chain whose send cache was empty (never solved) the correlator
    re-solved the plain ground state over the generalized one: |<wg|wf0>|^2
    = 0.60, <Sz0> -0.48 -> 0, while a chain solved first read wg. The
    metric A = 1 + 0.8*Sz0 is positive definite and breaks the Sz0 -> -Sz0
    symmetry, so the two states differ at order one."""
    np.random.seed(3)
    sc = _chain(version, n=6, maxm=20, nsweeps=10)
    if solved_first: _quiet(sc.gs_energy)
    lam = _quiet(lambda: sc.gs_energy_generalized(1 + 0.8*sc.Sz[0]))
    wg = sc.wf0.copy()
    assert groundstate.hamiltonian_on_session(sc)
    _quiet(lambda: sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[1]),
                                               es=np.linspace(-1, 4, 11)))
    assert _fidelity(sc.wf0, wg) == pytest.approx(1.0, abs=1e-10)
    assert sc.e0 == pytest.approx(lam, abs=1e-12)
    assert np.real(sc.vev(sc.Sz[0])) == pytest.approx(
        np.real(_expect(wg, sc.Sz[0])), abs=1e-10)
    assert abs(np.real(_expect(wg, sc.Sz[0]))) > 0.1


@pytest.mark.parametrize("version", BACKENDS)
@pytest.mark.parametrize("solved_first", [False, True], ids=["fresh", "solved"])
def test_nh_correlator_after_gs_energy_generalized_reads_the_generalized_state(
        version, solved_first):
    """The non-Hermitian route: neither NH solve put H on the session, so
    the next NH-KPM re-solved plain NH-DMRG over the generalized state on
    a fresh chain and on a solved one alike (|<wg|wf0>|^2 = 0.6877, e0 from
    lambda to -1.596396, <Sz0> -0.4529 -> 0)."""
    np.random.seed(9)
    sc = _nh_chain(version)
    if solved_first: _quiet(sc.gs_energy)
    lam = _quiet(lambda: sc.gs_energy_generalized(1 + 0.8*sc.Sz[0]))
    wg = sc.wf0.copy()
    assert groundstate.hamiltonian_on_session(sc)
    _quiet(lambda: sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]), **NHKPM))
    assert _fidelity(sc.wf0, wg) == pytest.approx(1.0, abs=1e-10)
    assert sc.e0 == pytest.approx(lam, abs=1e-12)
    assert sc.vev(sc.Sz[0]) == pytest.approx(_expect(wg, sc.Sz[0]), abs=1e-10)
    assert abs(np.real(_expect(wg, sc.Sz[0]))) > 0.1


@pytest.mark.parametrize("version", BACKENDS)
def test_first_nh_correlator_does_not_resolve_a_solved_chain(version,
                                                             monkeypatch):
    """After a plain NH-DMRG solve the first correlator solved NH-DMRG a
    second time, since the solve had not put H on the session; it now
    reads the solve's own pair."""
    from dmrgpy import nhdmrg
    calls = []
    real = nhdmrg.nhdmrg
    def counted(*args, **kwargs):
        calls.append(1)
        return real(*args, **kwargs)
    monkeypatch.setattr(nhdmrg, "nhdmrg", counted)
    np.random.seed(9)
    sc = _nh_chain(version)
    _quiet(sc.gs_energy)
    assert len(calls) == 1
    wf, left = sc.wf0, sc.nh_left_wf
    _quiet(lambda: sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]), **NHKPM))
    assert len(calls) == 1
    assert sc.wf0 is wf and sc.nh_left_wf is left


# --------------------------------------------------------- thermal-bypass

@pytest.fixture(scope="module")
def thermal3():
    """Exact Boltzmann <Sz0 Sz1> of the 3-site Heisenberg chain at T=1"""
    ref = spinchain.Spin_Chain(["S=1/2"]*3)
    ref.set_hamiltonian(_heis(ref))
    ed = ref.get_ED_obj()
    H = np.asarray(ed.get_hamiltonian().todense())
    ZZ = np.asarray(ed.MO2matrix(ref.Sz[0]*ref.Sz[1]).todense())
    w, U = np.linalg.eigh(H)
    p = np.exp(-(w - w[0])); p = p/p.sum()
    return float(np.real(np.sum(p*np.diag(U.conj().T @ ZZ @ U))))


@pytest.mark.parametrize("version", BACKENDS + [_backend(2), pytest.param("ED", id="ED")])
def test_every_reader_of_the_thermal_chain_measures_the_annealed_state(
        thermal3, version):
    """MBChain.vev() and the KPM sum rule on mode="ED" read the singlet
    state (0.0 against -0.0694), MBChain.gs_energy() on a DMRG backend
    returned the singlet Hamiltonian's -2.25 against <wf|H|wf> = -0.4164,
    and a KPM correlator re-solved the plain ground state (sum rule -0.1664,
    vev afterwards -0.1667). Now vev, the correlator's sum rule and, on the
    DMRG backends, gs_energy() all read the annealed state, whose <Sz0 Sz1>
    is the exact Boltzmann -0.0714 up to anneal()'s Euler steps. On
    mode="ED" gs_energy() is the lowest eigenvalue, which the 2026-09-24c
    record leaves as it is, so it is not pinned here."""
    np.random.seed(4)
    kw = {} if version == "ED" else dict(itensor_version=version)
    tc = thermal.Thermal_Spin_Chain(["S=1/2"]*3, T=1.0, **kw)
    ht = 0
    for i in range(2):
        ht = ht + tc.Sx[i]*tc.Sx[i+1] + tc.Sy[i]*tc.Sy[i+1] + tc.Sz[i]*tc.Sz[i+1]
    tc.set_hamiltonian(ht)
    if version == "ED": tc.mode = "ED"
    mb = tc.MBChain
    if version != "ED": mb.maxm, mb.nsweeps = 30, 10
    wf = _quiet(tc.get_gs)
    zz = tc.Sz[0]*tc.Sz[1]
    own = np.real(_expect(wf, zz))
    assert own == pytest.approx(thermal3, abs=5e-3)
    assert np.real(_quiet(lambda: mb.vev(zz))) == pytest.approx(own, abs=1e-10)
    if version != "ED":
        assert _quiet(mb.gs_energy) == pytest.approx(_energy(mb, wf), abs=1e-10)
    kwm = dict(mode="ED") if version == "ED" else {}
    es = np.linspace(-6.0, 6.0, 601)
    es, d = _quiet(lambda: mb.get_dynamical_correlator(
        name=(tc.Sz[0], tc.Sz[1]), submode="KPM", es=es, delta=0.2, **kwm))
    assert float(np.real(np.trapezoid(d, es))) == pytest.approx(own, abs=1e-3)
    if version != "ED":
        assert _fidelity(mb.wf0, wf) == pytest.approx(1.0, abs=1e-10)
    assert np.real(_quiet(lambda: mb.vev(zz))) == pytest.approx(own, abs=1e-10)


@pytest.mark.parametrize("version", BACKENDS)
def test_zero_temperature_thermal_chain_is_the_solved_ground_state(version):
    """The T<=1e-5 branch no longer writes a copy back over MBChain's own
    solved state; the state returned is that ground state, normalized."""
    np.random.seed(5)
    tc = thermal.Thermal_Spin_Chain(["S=1/2"]*3, T=0.0, itensor_version=version)
    ht = tc.Sx[0]*tc.Sx[1] + tc.Sy[0]*tc.Sy[1] + tc.Sz[0]*tc.Sz[1]
    tc.set_hamiltonian(ht)
    tc.MBChain.maxm, tc.MBChain.nsweeps = 20, 10
    wf = _quiet(tc.get_gs)
    assert np.real(wf.dot(wf)) == pytest.approx(1.0, abs=1e-12)
    assert _energy(tc.MBChain, wf) == pytest.approx(-0.75, abs=1e-8)
    assert _fidelity(tc.MBChain.get_gs(), wf) == pytest.approx(1.0, abs=1e-10)


# ----------------------------------------------------------- nh-injected

@pytest.mark.parametrize("version", BACKENDS)
def test_nh_kpm_refuses_a_right_state_with_no_left_partner(version):
    """After set_gs(x) the NH-KPM paired x with the left state of the last
    NH-DMRG solve (<psil|x> = 0.93 on this chain, and a spectrum with it),
    and on a chain that never ran NH-DMRG it raised AttributeError."""
    np.random.seed(5)
    sc = _nh_chain(version)
    _quiet(sc.gs_energy)
    psir = sc.wf0.copy()
    sc.set_gs(_normalized(psir + 0.5*(sc.Sx[2]*psir)))
    with pytest.raises(RuntimeError, match="left eigenvector"):
        _quiet(lambda: sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]), **NHKPM))
    fresh = _nh_chain(version)
    fresh.set_gs(_normalized(fresh.random_state()))
    with pytest.raises(RuntimeError, match="left eigenvector"):
        _quiet(lambda: fresh.get_dynamical_correlator(name=(fresh.Sz[0], fresh.Sz[0]),
                                                      **NHKPM))


@pytest.mark.parametrize("version", BACKENDS)
def test_nh_kpm_runs_on_the_pair_of_a_solve(version):
    """The pairing check passes on the chain's own NH-DMRG pair, before and
    after a set state is replaced by a fresh solve."""
    np.random.seed(6)
    sc = _nh_chain(version)
    x, y = _quiet(lambda: sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]), **NHKPM))
    assert np.all(np.isfinite(y))
    sc.set_gs(_normalized(sc.random_state()))
    sc.restart()
    x, y2 = _quiet(lambda: sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]), **NHKPM))
    assert np.max(np.abs(y2 - y)) < 1e-3*np.max(np.abs(y))


@pytest.mark.parametrize("version", BACKENDS + [_backend(2)])
def test_nh_gs_energy_wf0_is_taken_or_refused(version):
    """gs_energy(wf0=x) on a non-Hermitian chain returned the NH-DMRG energy
    with |<x|wf0>|^2 = 0.02 to 0.11: NH-DMRG takes no start state, so x is
    refused, and taken as it is with reconverge=False."""
    np.random.seed(7)
    sc = _nh_chain(version)
    x = _normalized(sc.random_state())
    with pytest.raises(TypeError, match="no start state"):
        sc.gs_energy(wf0=x)
    e = sc.gs_energy(wf0=x, reconverge=False)
    assert e == pytest.approx(_expect(x, sc.hamiltonian), abs=1e-10)
    assert _fidelity(sc.wf0, x) == pytest.approx(1.0, abs=1e-12)


# -------------------------------------------------------- maxde-per-site

def test_maxde_is_documented_per_site():
    for f in (Many_Body_Chain.gs_energy, groundstate.gs_energy_single):
        doc = f.__doc__ or ""
        assert "maxde" in doc and re.search(r"per.site", doc, re.I)
    assert re.search(r"divided by the number\s+of sites",
                     Many_Body_Chain.gs_energy_fluctuation.__doc__)


@pytest.mark.parametrize("version", BACKENDS)
def test_maxde_compares_the_fluctuation_per_site(version):
    """The loop's first reading is gs_energy_fluctuation()/ns of the first
    solve, and it stops once gs_energy_fluctuation()/ns is below maxde."""
    np.random.seed(8)
    sc = _chain(version, n=8, maxm=3, nsweeps=10)
    _quiet(sc.gs_energy)
    total = sc.gs_energy_fluctuation()
    maxde = 0.5*total/sc.ns # between the per-site and the total value
    # the next read takes the session's own solved state and its cached
    # energy, so the loop starts from the state just measured
    sc.computed_gs = False
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        sc.gs_energy(maxde=maxde)
    reads = [float(v) for v in re.findall(
        r"Energy fluctuation per site =\s+(\S+)", buf.getvalue())]
    assert reads and reads[0] == pytest.approx(total/sc.ns, rel=1e-3)
    assert sc.gs_energy_fluctuation()/sc.ns < maxde
