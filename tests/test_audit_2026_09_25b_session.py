"""Regression tests for the `session` cluster of the 2026-09-25b fix pass
(docs/audit_2026_09_25b_hole_hunt.md, findings 1, 8, 9 and 10, and the
julia_live solver-key lead of docs/audit_2026_09_25_open_items.md).

  finding 1. A state set with a norm other than one was a ray to some
      readers and a vector to others: after set_gs(2*s) gs_energy() and the
      session vev() divided by <x|x>, while KPM, CVM, TD, TDZ, ROOTN and
      evolve_and_measure() did, and on mode="ED" and julia_live vev() did
      not either, so the KPM sum rule came back 4 times the chain's own
      vev() and gs_energy_fluctuation() of an exact eigenstate was 6|E0|.
      A state is now normalized once, where it becomes the chain's
      (groundstate.unit_copy()), and the zero state is refused.
  finding 8. After gs_energy_generalized(A) e0 was lambda, from which KPM,
      CVM, ROOTN and TDZ measured the generalized state wg while TD and EX
      measured it from <wg|H|wg>: for A = 2*Id, where wg is the plain
      ground state, the KPM put its lines at -0.15 and +0.56 instead of
      +0.66 and +1.365. e0 is now wg's own energy (the biorthogonal one on
      the non-Hermitian route) and lambda is returned and kept as
      lam_generalized.
  finding 9. On a non-Hermitian chain gs_energy(H=H2) stored H2's pair as
      the chain's state, and the next NH-KPM of the chain's own H1 read it,
      0.587 of the peak off. H= is now refused, as on the Hermitian route.
  finding 10. The julia_live KPM window took its lower edge from e0 unless
      the state was marked supplied, so after gs_energy_generalized() it
      was lambda, and A = 2*Id or 1.5+0.4*Sz0 raised "KPM moments
      diverging". It now takes e0 only when a plain solve marked it as H's
      lower edge (groundstate.mark_lower_edge()).
  julia solver key. julia_live recorded no solver key, so a maxm ramp on
      one chain returned the first energy every time.

Anchors: ray invariance against the unit vector itself, exact expectation
values, ED (dense eigh and eig of the generalized problems), and the plain
ground state's own lines for the scalar metric.
"""

import io
import contextlib
import warnings

import numpy as np
import pytest
import scipy.linalg as sla

from dmrgpy import cppext, spinchain, timedependent

from _helpers import julia_available


def _backend(version):
    marks = ()
    if version in (2, 3):
        marks = pytest.mark.skipif(not cppext.available(version),
            reason="needs the compiled itensor_version=%d extension" % version)
    ident = "v%d" % version if version in (2, 3) else str(version)
    return pytest.param(version, id=ident, marks=marks)


SESSION = [_backend("python"), _backend(3)]
DMRG = SESSION + [_backend(2)]


def _require_julia():
    """Deferred to the test body so that collecting this file (and
    -k "not julia_live") never boots a Julia session."""
    ok, reason = julia_available()
    if not ok:
        pytest.skip("requires a working juliacall/Julia toolchain: %s" % reason)


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


def _chain(version, n=4, B=0.0, maxm=16, nsweeps=12):
    kw = {} if version == "ED" else dict(itensor_version=version)
    sc = spinchain.Spin_Chain(["S=1/2"]*n, **kw)
    sc.set_hamiltonian(_heis(sc, B))
    sc.maxm, sc.nsweeps = maxm, nsweeps
    return sc


def _expect(wf, op, left=None):
    bra = wf if left is None else left
    return complex(bra.dot(op*wf)/bra.dot(wf))


def _fidelity(a, b):
    return abs(a.dot(b))**2/abs(a.dot(a)*b.dot(b))


def _normalized(wf):
    return wf*(1/np.sqrt(wf.dot(wf).real))


def _peaks(grid, y, frac=0.05):
    y = np.real(y); m = np.max(y)
    return [float(grid[k]) for k in range(1, len(y)-1)
            if y[k] >= y[k-1] and y[k] > y[k+1] and y[k] > frac*m]


def _lines_match(found, exact, tol):
    return len(found) == len(exact) and all(
        abs(a-b) <= tol for a, b in zip(found, exact))


# ---------------------------------------------------------------- finding 1

ES1 = np.linspace(-1.0, 4.0, 501)


def _solved_unit_state(version):
    np.random.seed(1)
    sc = _chain(version)
    if version == "ED":
        s = _quiet(lambda: sc.get_gs(mode="ED"))
    else:
        _quiet(sc.gs_energy)
        s = _normalized(sc.get_gs().copy())
    return sc, s


def _ray_readers(sc, ed):
    kw = dict(mode="ED") if ed else {}
    _, y = _quiet(lambda: sc.get_dynamical_correlator(
        name=(sc.Sz[0], sc.Sz[0]), submode="KPM", es=ES1, delta=0.1, **kw))
    w = _quiet(lambda: sc.get_gs(**kw))
    return dict(
        vev=np.real(_quiet(lambda: sc.vev(sc.Sz[0]*sc.Sz[1], **kw))),
        fluct=np.real(_quiet(lambda: sc.gs_energy_fluctuation(**kw))),
        kpm=float(np.real(np.trapezoid(y, ES1))),
        norm2=float(np.real(w.dot(w))))


def _set(sc, route, x):
    if route == "set_gs": sc.set_gs(x)
    elif route == "set_initial_wf": sc.set_initial_wf(x)
    elif route == "gs_energy(wf0=)":
        _quiet(lambda: sc.gs_energy(wf0=x, reconverge=False))
    else: _quiet(lambda: sc.get_gs(wf0=x, reconverge=False))


@pytest.mark.parametrize("route", ["set_gs", "set_initial_wf",
                                   "gs_energy(wf0=)", "get_gs(wf0=)"])
@pytest.mark.parametrize("version", DMRG)
def test_a_set_state_is_the_ray_it_names(version, route):
    """set_gs(2*s) of the solved ground state s of a 4-site Heisenberg
    chain, and the three other routes that make a state the chain's
    unswept: every reader gives the number of s itself. The KPM sum rule
    was 0.999995 against <s|Sz0 Sz0|s> = 0.25 on "python", v3 and v2, and
    get_gs() handed back <w|w> = 4, while vev() and the fluctuation, which
    the session normalizes, were right already."""
    sc, s = _solved_unit_state(version)
    ref = dict(vev=np.real(_expect(s, sc.Sz[0]*sc.Sz[1])), fluct=0.0,
               kpm=0.25, norm2=1.0)
    for c in (2.0, 0.5):
        _set(sc, route, s*c)
        got = _ray_readers(sc, ed=False)
        for k in ref:
            assert got[k] == pytest.approx(ref[k], abs=1e-4), (c, k)


def test_a_set_state_is_the_ray_it_names_on_ed():
    """mode="ED": set_gs(2*s) stored 2*s as given, so every reader, vev()
    included, came back 4 times s's (vev(Sz0 Sz1) -0.910684 against
    -0.227671, the fluctuation 9.696152 = 6|E0| for an exact eigenstate),
    and so did the correlator's own explicit wf0=."""
    sc, s = _solved_unit_state("ED")
    ref = dict(vev=np.real(_expect(s, sc.Sz[0]*sc.Sz[1])), fluct=0.0,
               kpm=0.25, norm2=1.0)
    for c in (2.0, 0.5):
        sc.set_gs(s*c)
        got = _ray_readers(sc, ed=True)
        for k in ref:
            assert got[k] == pytest.approx(ref[k], abs=1e-4), (c, k)
    grid = np.linspace(-1.0, 4.0, 201)
    out = []
    for c in (1.0, 2.0):
        _, y = _quiet(lambda: sc.get_dynamical_correlator(mode="ED",
            name=(sc.Sz[0], sc.Sz[0]), submode="CVM", es=grid, delta=0.1,
            wf0=s*c))
        out.append(np.asarray(y))
    assert np.max(np.abs(out[1] - out[0])) < 1e-10*np.max(np.abs(out[0]))


@pytest.mark.parametrize("version", DMRG)
def test_evolve_and_measure_reads_the_ray_but_not_an_explicit_start(version):
    """evolve_and_measure() without wf= reads self.wf0 before any
    ground-state read, which is why the normalization sits where the state
    is marked and not only where it is taken: straight after set_gs(2*s)
    its t=0 value was 4*<s|Sz0 Sz1|s> (-0.910684). An explicit wf= is left
    as given, <psi(t)|O|psi(t)> of the raw vector, which is that route's
    own contract (2026-08 audit, finding 9)."""
    sc, s = _solved_unit_state(version)
    op = sc.Sz[0]*sc.Sz[1]
    ref = np.real(_expect(s, op))
    sc.set_gs(s*2.0)
    _, cs = _quiet(lambda: timedependent.evolve_and_measure_dmrg(
        sc, operator=op, nt=2, dt=0.01))
    assert np.real(cs[0]) == pytest.approx(ref, abs=1e-6)
    _, cs = _quiet(lambda: timedependent.evolve_and_measure_dmrg(
        sc, operator=op, nt=2, dt=0.01, wf=s*2.0))
    assert np.real(cs[0]) == pytest.approx(4*ref, abs=1e-6)


@pytest.mark.parametrize("version", SESSION + ["ED"])
def test_the_zero_state_is_refused(version):
    """normalize() returns None below its floor, so the fix must refuse a
    state with no direction rather than store None."""
    sc, s = _solved_unit_state(version)
    with pytest.raises(ValueError, match="no direction"):
        sc.set_gs(s*0.0)
    if version != "ED":
        with pytest.raises(ValueError, match="no direction"):
            sc.gs_energy(wf0=s*0.0, reconverge=False)
        with pytest.raises(ValueError, match="no direction"):
            sc.set_initial_wf_guess(s*0.0)


# ----------------------------------------------------------- findings 8, 10

ES8 = np.linspace(-1.5, 3.5, 1001)
METRICS = {"2Id": lambda sc: 8*sc.Sz[0]*sc.Sz[0],
           "1.5+0.4Sz0": lambda sc: 1.5 + 0.4*sc.Sz[0],
           "1+0.8Sz0": lambda sc: 1 + 0.8*sc.Sz[0]}


def _exact_generalized(which):
    """lam, E_wg = <v|H|v> and the (Sz0,Sz0) lines E_n - E_wg of weight
    above 0.005 (the next one down is 0.0036), from a dense eigh(H, A) of
    the 4-site chain in Bz=0.3."""
    ref = _chain("ED", B=0.3)
    ed = ref.get_ED_obj()
    H = np.asarray(ed.get_hamiltonian().todense())
    A = np.asarray(ed.MO2matrix(METRICS[which](ref)).todense())
    Z0 = np.asarray(ed.MO2matrix(ref.Sz[0]).todense())
    lams, vecs = sla.eigh(H, A)
    v = vecs[:, 0]/np.linalg.norm(vecs[:, 0])
    ewg = float(np.real(np.vdot(v, H@v)))
    E, U = np.linalg.eigh(H)
    w = np.abs(U.conj().T @ (Z0@v))**2
    return lams[0], ewg, E[0], sorted(E[k]-ewg for k in range(len(E)) if w[k] > 0.005)


def _generalized_chain(version, which):
    np.random.seed(2)
    sc = _chain(version, B=0.3)
    lam = _quiet(lambda: sc.gs_energy_generalized(METRICS[which](sc)))
    return sc, lam, sc.wf0.copy()


@pytest.mark.parametrize("which", list(METRICS))
@pytest.mark.parametrize("version", SESSION)
def test_generalized_state_is_measured_from_its_own_energy(version, which):
    """After gs_energy_generalized(A) the KPM, which measures from e0, put
    wg's lines lambda - <wg|H|wg> away from where TD and EX put them: on
    1+0.8*Sz0 at [0.545, 1.205, 1.915] against the exact E_n - E_wg of
    [-0.230, 0.429, 1.136]. Now it lands on E_n - E_wg, gs_energy() and e0
    are E_wg (gs_energy() returned lambda), and lambda is still the return
    value and lam_generalized. A = 1.5+0.4*Sz0 and 2*Id put lambda above
    E0 with weight at E0, the metrics that made julia_live's window raise
    (finding 10): the session backends' own window solves for H's lower
    edge, and holds."""
    lam_ex, ewg, _, lines = _exact_generalized(which)
    sc, lam, wg = _generalized_chain(version, which)
    assert lam == pytest.approx(lam_ex, abs=1e-6)
    assert sc.lam_generalized == lam
    assert sc.e0 == pytest.approx(ewg, abs=1e-6)
    assert _quiet(sc.gs_energy) == pytest.approx(ewg, abs=1e-6)
    _, y = _quiet(lambda: sc.get_dynamical_correlator(
        name=(sc.Sz[0], sc.Sz[0]), submode="KPM", es=ES8, delta=0.05))
    assert _lines_match(_peaks(ES8, y), lines, 0.0051), _peaks(ES8, y)
    assert np.real(np.trapezoid(y, ES8)) == pytest.approx(0.25, abs=1e-3)
    assert _fidelity(sc.wf0, wg) == pytest.approx(1.0, abs=1e-10)


@pytest.mark.parametrize("version", SESSION)
def test_scalar_metric_puts_every_submode_on_the_plain_ground_state_lines(version):
    """A = 2*Id, where the generalized state is exactly the plain ground
    state and lambda = E0/2: ED puts the (Sz0,Sz0) lines at E_n - E0 =
    +0.659 and +1.366, where KPM, CVM, ROOTN and TDZ put them at E_n -
    lambda = -0.149 and +0.558 (a line at negative frequency in a
    ground-state autocorrelator) and TD and EX at the right place. CVM and
    ROOTN cost a solve per frequency, so they are probed at the four
    positions only."""
    _, ewg, e0_plain, lines = _exact_generalized("2Id")
    assert ewg == pytest.approx(e0_plain, abs=1e-10)
    assert np.round(lines, 3).tolist() == [0.659, 1.366]
    sc, lam, wg = _generalized_chain(version, "2Id")
    assert lam == pytest.approx(e0_plain/2, abs=1e-6)
    name = (sc.Sz[0], sc.Sz[0])
    for sub in ("KPM", "TD", "TDZ", "EX"):
        kw = dict(nex=16) if sub == "EX" else {}
        _, y = _quiet(lambda: sc.get_dynamical_correlator(
            name=name, submode=sub, es=ES8, delta=0.05, **kw))
        assert _lines_match(_peaks(ES8, y), lines, 0.0051), (sub, _peaks(ES8, y))
    wrong = [e0_plain - lam + d for d in lines] # the lines from lambda
    probe = np.array([wrong[0], wrong[1], lines[0], lines[1]])
    for sub in ("CVM", "ROOTN"):
        _, y = _quiet(lambda: sc.get_dynamical_correlator(
            name=name, submode=sub, es=probe, delta=0.05))
        y = np.real(np.asarray(y))
        # the exact lines against the lambda ones: a Lorentzian of width
        # 0.05 is 5 times lower 0.1 off its line and 100 times 0.8 off
        assert y[2] > 3*y[1] and y[3] > 20*y[0], (sub, y)
    assert _fidelity(sc.wf0, wg) == pytest.approx(1.0, abs=1e-10)


@pytest.mark.parametrize("version", SESSION)
def test_generalized_state_reads_as_the_same_state_set_by_hand(version):
    """The same wg two ways onto one chain, as gs_energy_generalized()'s
    result and handed back with set_gs(wg): the KPM curve and gs_energy()
    were 0.777 apart (lambda against <wg|H|wg>) and are now the same."""
    np.random.seed(2)
    sc = _chain(version, B=0.3)
    _quiet(lambda: sc.gs_energy_generalized(1 + 0.8*sc.Sz[0]))
    wg = sc.wf0.copy()
    out = []
    for step in range(2):
        if step: sc.set_gs(wg)
        _, y = _quiet(lambda: sc.get_dynamical_correlator(
            name=(sc.Sz[0], sc.Sz[0]), submode="KPM", es=ES8, delta=0.05))
        out.append((np.asarray(y), _quiet(sc.gs_energy)))
    assert out[1][1] == pytest.approx(out[0][1], abs=1e-10)
    assert np.max(np.abs(out[1][0] - out[0][0])) < 1e-6*np.max(np.abs(out[0][0]))


def _nh_chain(version):
    sc = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version=version)
    sc.set_hamiltonian(_heis(sc) + 0.3j*sc.Sz[0] + 0.2*sc.Sx[1])
    sc.maxm, sc.nsweeps = 20, 10
    return sc


@pytest.mark.parametrize("version", SESSION)
def test_nh_generalized_energy_is_the_pairs_own(version):
    """The non-Hermitian route stored e0 = lambda too. It is now the pair's
    own biorthogonal energy <psil|H|psir>/<psil|psir>, anchored on the
    exact left and right generalized eigenvectors of a dense eig(H, A)."""
    ref = _nh_chain("python")
    ed = ref.get_ED_obj()
    H = np.asarray(ed.get_hamiltonian().todense())
    A = np.asarray(ed.MO2matrix(1 + 0.8*ref.Sz[0]).todense())
    w, vl, vr = sla.eig(H, A, left=True, right=True)
    k = np.argmin(w.real)
    u, v = vl[:, k], vr[:, k]
    e_exact = np.vdot(u, H@v)/np.vdot(u, v)
    np.random.seed(9)
    sc = _nh_chain(version)
    lam = _quiet(lambda: sc.gs_energy_generalized(1 + 0.8*sc.Sz[0]))
    assert lam == pytest.approx(w[k], abs=1e-5)
    assert sc.lam_generalized == lam
    assert abs(e_exact - w[k]) > 0.1 # the two are different numbers here
    assert sc.e0 == pytest.approx(e_exact, abs=1e-5)
    assert sc.e0 == pytest.approx(_expect(sc.wf0, sc.hamiltonian,
                                          left=sc.nh_left_wf), abs=1e-10)
    assert _quiet(sc.gs_energy) == sc.e0


# ---------------------------------------------------------------- finding 9

@pytest.mark.parametrize("version", SESSION)
def test_nh_gs_energy_refuses_another_operator(version):
    """gs_energy(H=H2) on the chain of H1 stored H2's pair, and after one
    solve and set_initial_wf(None) the next NH-KPM of H1 read H2's pair and
    energy (-1.836506+0.072051j against -1.596396), 0.587 of the peak off.
    H= is refused now, as on the Hermitian route; nhdmrg(H=...) still
    returns the other operator's pair without storing it."""
    def H2(sc): return _heis(sc) + 0.3j*sc.Sz[0] + 1.0*sc.Sz[3]
    KW = dict(submode="KPM", es=np.linspace(0.0, 3.0, 13), delta=0.3,
              E_max=10.0, n=60)
    np.random.seed(3)
    sc = _nh_chain(version)
    _quiet(sc.gs_energy)
    sc.set_initial_wf(None)
    with pytest.raises(TypeError, match="nhdmrg"):
        _quiet(lambda: sc.gs_energy(H=H2(sc)))
    _, y = _quiet(lambda: sc.get_dynamical_correlator(name=(sc.Sz[3], sc.Sz[3]), **KW))
    e = _quiet(sc.gs_energy)
    sc.restart()
    _, yref = _quiet(lambda: sc.get_dynamical_correlator(name=(sc.Sz[3], sc.Sz[3]), **KW))
    assert e == pytest.approx(-1.596396, abs=1e-5)
    assert np.max(np.abs(np.asarray(y) - yref)) < 1e-6*np.max(np.abs(yref))
    fresh = _nh_chain(version)
    with pytest.raises(TypeError, match="nhdmrg"):
        _quiet(lambda: fresh.gs_energy(H=H2(fresh)))
    e2, _, _ = _quiet(lambda: fresh.nhdmrg(H=H2(fresh)))
    assert e2 == pytest.approx(-1.836506+0.072051j, abs=1e-5)
    assert not fresh.computed_gs


# --------------------------------------------------------------- julia_live

@pytest.mark.parametrize("version", [pytest.param("julia_live", id="julia_live")])
def test_julia_set_state_is_the_ray_it_names(version):
    """julia_live's vev(), fluctuation, KPM and excited states read the set
    vector as given: after set_gs(2*s) vev(Sz0 Sz1) was -0.910684, the
    fluctuation 9.696152, the KPM sum rule 0.999995 and the excited-state
    search reported the set state at 4*E0 = -6.464102. Normalizing at
    injection covers them without a change of their own."""
    _require_julia()
    sc, s = _solved_unit_state("julia_live")
    e0 = _quiet(sc.gs_energy)
    ref = np.real(_expect(s, sc.Sz[0]*sc.Sz[1]))
    sc.set_gs(s*2.0)
    got = _ray_readers(sc, ed=False)
    assert got["vev"] == pytest.approx(ref, abs=1e-6)
    assert got["fluct"] == pytest.approx(0.0, abs=1e-3)
    assert got["kpm"] == pytest.approx(0.25, abs=1e-4)
    assert got["norm2"] == pytest.approx(1.0, abs=1e-10)
    ee, _ = _quiet(lambda: sc.get_excited_states(n=2))
    assert np.real(ee[0]) == pytest.approx(np.real(e0), abs=1e-6)


@pytest.mark.parametrize("version", [pytest.param("julia_live", id="julia_live")])
def test_julia_generalized_kpm_measures_wg_on_the_band_edge_window(version):
    """julia_live's KPM window took its lower edge from e0 = lambda after a
    generalized solve, so A = 2*Id and A = 1.5+0.4*Sz0 raised "KPM moments
    diverging", and on 1+0.8*Sz0 (lambda below E0) the window was too wide
    and the line 0.156 of the peak off the session backends'. Now the edge
    is a solve of H's own, the origin is <wg|H|wg> (finding 8), and the
    curve is "python"'s on the same call."""
    _require_julia()
    for which in METRICS:
        lam_ex, ewg, _, lines = _exact_generalized(which)
        curves = []
        for v in ("julia_live", "python"):
            np.random.seed(2)
            sc = _chain(v, B=0.3)
            lam = _quiet(lambda: sc.gs_energy_generalized(METRICS[which](sc)))
            wg = sc.wf0.copy()
            assert lam == pytest.approx(lam_ex, abs=1e-6)
            assert sc.e0 == pytest.approx(ewg, abs=1e-6)
            _, y = _quiet(lambda: sc.get_dynamical_correlator(
                name=(sc.Sz[0], sc.Sz[0]), submode="KPM", es=ES8, delta=0.05))
            curves.append(np.asarray(y))
            assert _lines_match(_peaks(ES8, y), lines, 0.0051), (which, v)
            assert np.real(np.trapezoid(y, ES8)) == pytest.approx(0.25, abs=1e-3)
            assert _fidelity(sc.wf0, wg) == pytest.approx(1.0, abs=1e-8)
            assert sc.e0 == pytest.approx(ewg, abs=1e-6)
        assert np.max(np.abs(curves[0] - curves[1])) < 1e-2*np.max(np.abs(curves[1])), which


@pytest.mark.parametrize("version", [pytest.param("julia_live", id="julia_live")])
def test_julia_solver_key(version):
    """julia_live recorded no solver key, so a maxm ramp on one 8-site
    chain returned -3.194321 at maxm 2, 4 and 16. The key the band-edge
    solves of the KPM and of the excited-state search leave behind (they
    run at a clamped maxm/nsweeps) must not make a state look stale
    afterwards, or the next read replaces a set state with a solve."""
    _require_julia()
    np.random.seed(3)
    sk = _chain("julia_live", n=8, nsweeps=10)
    es = []
    for m in (2, 4, 16):
        sk.maxm = m
        es.append(np.real(_quiet(sk.gs_energy)))
    assert es[0] > es[1] + 1e-3 and es[1] > es[2] + 1e-3
    assert es[2] == pytest.approx(-3.374932, abs=1e-4)
    # at the default maxm=30/nsweeps=15, above the helpers' 20/5 clamp
    np.random.seed(4)
    sc = _chain("julia_live", n=4, maxm=30, nsweeps=15)
    _quiet(sc.gs_energy)
    g = sc.get_gs().copy()
    x = _normalized(g + 0.5*(sc.Sx[0]*g))
    ex = np.real(_expect(x, sc.hamiltonian))
    sc.set_gs(x)
    _quiet(lambda: sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]),
            submode="KPM", es=np.linspace(-1, 4, 51), delta=0.2))
    _quiet(lambda: sc.get_excited_states(n=2))
    assert _quiet(sc.gs_energy) == pytest.approx(ex, abs=1e-8)
    assert _fidelity(sc.get_gs(), x) == pytest.approx(1.0, abs=1e-10)
