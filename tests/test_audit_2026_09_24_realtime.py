"""Regressions for the real-time cluster of the 2026-09-24 hole hunt
(docs/audit_2026_09_24_hole_hunt.md, findings 7 to 10).

What each section pins:

- Finding 7. `timedependent._fourier_transform_correlator` evaluated the
  damped trapezoid sum with an FFT, on a grid of spacing 2*pi/(nt*dt),
  1.05*delta at the default time window, and interpolated linearly onto
  `es`. That interpolation was the whole 3.31e-02 residual TDZ and TD at
  predict=False carried on the audit's complex-hopping chain, which the
  records had put down to TDZ's complex-time contour. The sum is now
  evaluated at each requested frequency, so the tests hold every
  real-time route to its own accuracy, and the contour on its own, as a
  TDZ-against-TD difference that shares everything but the contour. The
  infinite-chain `sxt_to_skomega` stays on the FFT stage on purpose.
- Finding 8. `mode="ED"` TD re-solved the ground state with a randomly
  started eigsh inside every evolution, so on a degenerate ground state
  the two halves of a pair that takes two evolutions came from two
  different members of the manifold. Both now use the EDchain's cached
  state (or the `wf0=` the caller passes), and the shift of H is by -E_0.
- Finding 9. `submode="TDZ"` swallowed every unknown keyword in a
  `**kwargs` it never read, and had lost its check that the operators
  are symbolic.
- Finding 10. The lower-level `get_dynamical_correlator_MB` route dropped
  `i=`/`j=` for TD and TDZ and answered for sites 0 and 0.
"""
import numpy as np
import pytest

from dmrgpy import fermionchain, spinchain, timedependent
from dmrgpy.edtk import timedependent as tded
from dmrgpy.edtk.edchain import State
from scipy.sparse import linalg as slg
from scipy.sparse import identity
from dmrgpy.edtk.tdtk import evolve


def _density(D, M, es, delta):
    """sum_n M_n L_delta(w - D_n), the house convention"""
    es = np.asarray(es)
    return ((delta / np.pi) / ((es[:, None] - D[None, :]) ** 2 + delta ** 2)) @ M


# ------------------------------------------------------------- finding 7

ES = np.linspace(-1.0, 6.0, 60)
DELTA, DT = 0.4, 0.1


def _complex_hopping_chain(n=4, seed=3):
    """The audit's seeded 4-site complex-hopping chain, and the exact
    Lehmann data of (Cdag_0, C_2) from a hand-built Jordan-Wigner kron,
    with no dmrgpy code in the reference."""
    fc = fermionchain.Fermionic_Chain(n)
    rng = np.random.RandomState(seed)
    t = rng.random((n, n)) + 1j * rng.random((n, n))
    t = t + t.conj().T
    h = 0
    for i in range(n):
        for j in range(n):
            h = h + t[i, j] * fc.Cdag[i] * fc.C[j]
    for i in range(n - 1):
        h = h + 0.8 * fc.N[i] * fc.N[i + 1]
    fc.set_hamiltonian(h)
    fc.maxm, fc.nsweeps = 30, 12
    a = np.array([[0, 1], [0, 0]], dtype=complex)
    Z = np.diag([1.0, -1.0]).astype(complex)
    I2 = np.eye(2, dtype=complex)

    def c(i):
        out = np.array([[1.0 + 0j]])
        for k in range(n):
            out = np.kron(out, Z if k < i else (a if k == i else I2))
        return out
    C = [c(i) for i in range(n)]
    Cd = [x.conj().T for x in C]
    Nn = [Cd[i] @ C[i] for i in range(n)]
    H = sum(t[i, j] * Cd[i] @ C[j] for i in range(n) for j in range(n))
    H = H + sum(0.8 * Nn[i] @ Nn[i + 1] for i in range(n - 1))
    e, U = np.linalg.eigh(H)
    Uh = U.conj().T
    M = (Uh @ Cd[0] @ U)[0, :] * (Uh @ C[2] @ U)[:, 0]
    return fc, e - e[0], M


@pytest.fixture(scope="module")
def audit_chain():
    return _complex_hopping_chain()


def _real_time(fc, mode, submode, **extra):
    _x, y = fc.get_dynamical_correlator(mode=mode, submode=submode,
                                        name=[fc.Cdag[0], fc.C[2]], es=ES,
                                        delta=DELTA, dt=DT, **extra)
    return np.asarray(y, dtype=np.complex128)


@pytest.mark.parametrize("mode,submode,extra,tol", [
    ("DMRG", "TDZ", {}, 1e-3),                  # was 3.31e-02, now 5.4e-04
    ("DMRG", "TD", {"predict": False}, 1e-3),   # was 3.31e-02, now 5.4e-04
    ("ED", "TD", {"predict": False}, 1e-3),     # was 3.31e-02, now 5.4e-04
    ("DMRG", "TD", {}, 1e-4),                   # was 2.57e-04, now 3.2e-05
    ("ED", "TD", {}, 1e-4),                     # was 2.57e-04, now 3.1e-05
])
def test_real_time_routes_are_held_to_their_own_accuracy(audit_chain, mode,
                                                         submode, extra, tol):
    """Pointwise against the exact complex Lehmann density, max|Im M_n| =
    0.27 on this pair, exact peak 0.2313. The residual left is the one
    set by the finite time window, not by how many grid points fall
    inside one line width."""
    fc, D, M = audit_chain
    assert np.max(np.abs(M.imag)) > 0.1, "this pair must have complex weights"
    err = np.max(np.abs(_real_time(fc, mode, submode, **extra)
                        - _density(D, M, ES, DELTA)))
    assert err < tol, "%s/%s %s is off the exact density by %.3e" \
        % (mode, submode, extra, err)


def test_the_contour_error_on_its_own(audit_chain):
    """TDZ and TD at predict=False share everything downstream of the
    time series, so their difference is the complex-time contour plus
    its Taylor-in-alpha0 reconstruction alone: measured 9.9e-07 on this
    chain, against the 3.31e-02 the records had attributed to it."""
    fc, _D, _M = audit_chain
    y_tdz = _real_time(fc, "DMRG", "TDZ")
    y_td = _real_time(fc, "DMRG", "TD", predict=False)
    assert np.max(np.abs(y_tdz - y_td)) < 3e-6


def _lehmann_series(nt, dt, seed=5):
    rng = np.random.RandomState(seed)
    D = rng.uniform(0.0, 4.0, 6)
    M = rng.normal(size=6) + 1j * rng.normal(size=6)
    ts = dt * np.arange(nt)
    return ts, np.exp(-1j * np.outer(ts, D)) @ M


@pytest.mark.parametrize("nt", [150, 151])
def test_direct_evaluation_is_the_fft_on_its_own_grid(nt, monkeypatch):
    """At a frequency on the FFT grid the direct sum and the FFT are the
    same number, which is what pins the time origin and the trapezoid
    weights of the direct evaluation; and chunking over `es` does not
    change the result."""
    dt = 0.1
    ts, cs = _lehmann_series(nt, dt)
    es = 2 * np.pi / (nt * dt) * np.arange(-3, 20)
    _e, y_fft = timedependent._fourier_transform_correlator(
        ts, cs, dt, es=es, delta=DELTA, _evaluation="fft")
    _e, y_direct = timedependent._fourier_transform_correlator(
        ts, cs, dt, es=es, delta=DELTA)
    assert np.max(np.abs(y_direct - y_fft)) < 1e-12
    monkeypatch.setattr(timedependent, "_DIRECT_FT_MAX_ELEMENTS", 7)
    _e, y_chunked = timedependent._fourier_transform_correlator(
        ts, cs, dt, es=es, delta=DELTA)
    assert np.max(np.abs(y_chunked - y_direct)) < 1e-12


def test_sxt_to_skomega_stays_on_the_fft_stage():
    """The infinite-chain S(k,omega) reduction was deliberately left on
    the FFT-plus-interpolation stage, since at its own defaults the series
    is cut at e^-1 of its envelope and nobody has measured what a direct
    evaluation does there. Pinned bit for bit, and in a regime where the
    two stages do differ, so the pin means something."""
    rng = np.random.RandomState(11)
    nt, nx, dt, delta = 25, 7, 0.05, 0.3
    ts = dt * np.arange(nt)
    xs = np.arange(nx) - nx // 2
    S = rng.normal(size=(nt, nx)) + 1j * rng.normal(size=(nt, nx))
    ks = np.linspace(-np.pi, np.pi, 5)
    es = np.linspace(-1.0, 6.0, 50)
    _ks, _es, Skw = timedependent.sxt_to_skomega(ts, xs, S, dt, ks=ks, es=es,
                                                  delta=delta)
    for ik, k in enumerate(ks):
        Skt = S @ np.exp(-1j * k * xs)
        _e, g_fft = timedependent._fourier_transform_correlator(
            ts, Skt, dt, es=es, delta=delta, _evaluation="fft")
        _e, g_direct = timedependent._fourier_transform_correlator(
            ts, Skt, dt, es=es, delta=delta)
        assert np.array_equal(Skw[ik], g_fft)
        assert np.max(np.abs(g_direct - g_fft)) > 1e-3


# ------------------------------------------------------------- finding 8

ES8 = np.linspace(-1.0, 3.0, 41)
DELTA8, DT8 = 0.3, 0.05


def _heisenberg3(hx=0.0, hz0=0.0):
    """3-site S=1/2 Heisenberg chain: a Kramers doublet ground state at
    hx=hz0=0, non-degenerate once the fields are on"""
    n = 3
    sc = spinchain.Spin_Chain([2] * n)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] \
                + sc.Sz[i] * sc.Sz[i + 1]
    for i in range(n):
        h = h + hx * sc.Sx[i]
    sc.set_hamiltonian(h + hz0 * sc.Sz[0])
    Sp = [sc.Sx[i] + 1j * sc.Sy[i] for i in range(n)]
    Sm = [sc.Sx[i] - 1j * sc.Sy[i] for i in range(n)]
    return sc, Sp, Sm


def _ed_matrices(sc, A, B):
    ed = sc.get_ED_obj()
    H = np.array(ed.get_hamiltonian().todense())
    Am = np.array(ed.MO2matrix(A).todense())
    Bm = np.array(ed.MO2matrix(B).todense())
    return ed, H, Am, Bm


def _state_density(H, Am, Bm, g):
    """The density of the pair in the state g, measured from E_0"""
    e, U = np.linalg.eigh(H)
    M = (g.conj() @ Am @ U) * (U.conj().T @ Bm @ g)
    return _density(e - e[0], M, ES8, DELTA8)


def _ed_td(sc, A, B, **kwargs):
    return np.asarray(sc.get_dynamical_correlator(
        mode="ED", submode="TD", name=[A, B], es=ES8, delta=DELTA8, dt=DT8,
        **kwargs)[1])


def _evolution_DC_before_the_fix(self, h=None, name=None, nt=100, dt=0.01,
                                 **kwargs):
    """edtk/timedependent.evolution_DC as it was, kept verbatim as the
    reference for the non-degenerate case (a fresh randomly started
    eigsh per call, and a shift by its eigenvalue of -H)"""
    (A, B) = name[1], name[0]
    Hop = self.get_operator(h)
    Aop = self.get_operator(A)
    Bop = self.get_operator(B)
    ts = np.array([dt * ii for ii in range(nt)])
    e0, wf0 = slg.eigsh(-Hop, k=1, ncv=20, which="LA")
    wf0 = wf0.reshape(wf0.shape[0])
    wf = wf0.copy()
    wf = Aop @ wf
    wfc = np.conjugate(wf0)
    cs = []
    ht = Hop + e0[0] * identity(Hop.shape[0], dtype=np.complex128)
    for it in range(nt):
        c = wfc @ Bop @ wf
        cs.append(c)
        wf = evolve(wf, ht, t=dt, dt=dt)
        wf = wf.reshape((wf.shape[0]))
    cs = np.array(cs)
    return ts, cs


def test_ed_td_on_a_degenerate_ground_state_uses_the_cached_state():
    """(Sp_0, Sm_2) is not provably self-adjoint, so its density takes two
    evolutions, and on this chain it depends on which member of the
    doublet it is measured in. Two identical calls used to differ by
    2.0e-01 on a peak of 1.8e-01, each off the cached-state density by
    1.4e-01 to 3.4e-01; now they are identical and within 3.3e-06 of it."""
    sc, Sp, Sm = _heisenberg3()
    ed, H, Am, Bm = _ed_matrices(sc, Sp[0], Sm[2])
    e = np.linalg.eigvalsh(H)
    assert e[1] - e[0] < 1e-10, "the ground state must be degenerate"
    y1 = _ed_td(sc, Sp[0], Sm[2])
    y2 = _ed_td(sc, Sp[0], Sm[2])
    assert np.array_equal(y1, y2)
    ref = _state_density(H, Am, Bm, ed.get_gs_array())
    assert np.max(np.abs(y1 - ref)) < 1e-4


def test_ed_td_measures_the_state_it_is_given():
    """wf0= reaches the evolution now: TD in another member of the
    doublet gives that member's density, not the cached one's."""
    sc, Sp, Sm = _heisenberg3()
    ed, H, Am, Bm = _ed_matrices(sc, Sp[0], Sm[2])
    _y = _ed_td(sc, Sp[0], Sm[2])  # computes and caches the ground state
    g0 = ed.get_gs_array()
    ref0 = _state_density(H, Am, Bm, g0)
    e, U = np.linalg.eigh(H)
    u0, u1 = U[:, 0], U[:, 1]
    s = 1 / np.sqrt(2)
    poles = [u0, u1, s * (u0 + u1), s * (u0 - u1), s * (u0 + 1j * u1),
             s * (u0 - 1j * u1)]
    refs = [_state_density(H, Am, Bm, g) for g in poles]
    k = int(np.argmax([np.max(np.abs(r - ref0)) for r in refs]))
    assert np.max(np.abs(refs[k] - ref0)) > 1e-2, \
        "the density must depend on the member for this test to mean much"
    y = _ed_td(sc, Sp[0], Sm[2], wf0=State(poles[k], ed))
    assert np.max(np.abs(y - refs[k])) < 1e-4


def test_ed_td_on_a_non_degenerate_ground_state_is_unchanged(monkeypatch):
    """Where the ground state is unique the fix changes nothing but
    rounding: against the pre-fix construction, same Fourier stage. It is
    also where the sign of the shift shows, on a pair that is not
    self-adjoint: keeping H + E_0 with the cached +E_0 would put the
    result 1.3e-01 off the exact density."""
    sc, Sp, Sm = _heisenberg3(hx=0.2, hz0=0.1)
    ed, H, Am, Bm = _ed_matrices(sc, Sp[0], Sm[2])
    e = np.linalg.eigvalsh(H)
    assert e[1] - e[0] > 1e-2, "the ground state must be non-degenerate"
    y_new = _ed_td(sc, Sp[0], Sm[2])
    monkeypatch.setattr(tded, "evolution_DC", _evolution_DC_before_the_fix)
    y_old = _ed_td(sc, Sp[0], Sm[2])
    monkeypatch.undo()
    assert np.max(np.abs(y_new - y_old)) < 1e-9
    ref = _state_density(H, Am, Bm, ed.get_gs_array())
    assert np.max(np.abs(y_new - ref)) < 1e-4


def test_evolution_aba_on_ed_uses_the_cached_state():
    """evolution_ABC, which backs evolution_ABA/evolve_and_measure on
    mode="ED", re-solved the ground state the same way when no wf= was
    given. It now uses the cached one: deterministic, and the same array
    as passing that state explicitly, while another member of the
    doublet gives a different one."""
    sc, Sp, Sm = _heisenberg3()
    ed = sc.get_ED_obj()
    kw = dict(A=Sp[0], B=sc.Sz[2], mode="ED", nt=40, dt=0.1)
    _t, c1 = timedependent.evolution_ABA(sc, **kw)
    _t, c2 = timedependent.evolution_ABA(sc, **kw)
    g0 = ed.get_gs_array()
    _t, c0 = timedependent.evolution_ABA(sc, wf=State(g0, ed), **kw)
    assert np.array_equal(c1, c2)
    assert np.array_equal(c1, c0)
    H = np.array(ed.get_hamiltonian().todense())
    _e, U = np.linalg.eigh(H)
    others = [U[:, 0], U[:, 1]]
    diffs = []
    for g in others:
        g = g - np.vdot(g0, g) * g0  # the member orthogonal to the cached one
        if np.linalg.norm(g) < 1e-6:
            continue
        g = g / np.linalg.norm(g)
        _t, cg = timedependent.evolution_ABA(sc, wf=State(g, ed), **kw)
        diffs.append(np.max(np.abs(np.asarray(cg) - np.asarray(c0))))
    assert max(diffs) > 1e-3


# ---------------------------------------------------- findings 9 and 10

@pytest.fixture(scope="module")
def heisenberg4():
    n = 4
    sc = spinchain.Spin_Chain([2] * n)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] \
                + sc.Sz[i] * sc.Sz[i + 1]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 20, 8
    return sc


BASE = dict(es=np.linspace(0.0, 3.0, 13), delta=0.4, dt=0.1)


@pytest.mark.parametrize("bad", [{"foo": 1}, {"alpha": 0.3}, {"nmax": 0}])
def test_tdz_rejects_an_unknown_keyword(heisenberg4, bad):
    """alpha= for alpha0= and nmax= for n_max= used to return the
    default-parameter spectrum bit for bit, where TD raised"""
    sc = heisenberg4
    with pytest.raises(TypeError) as excinfo:
        sc.get_dynamical_correlator(mode="DMRG", submode="TDZ",
                                    name=[sc.Sz[0], sc.Sz[0]], **BASE, **bad)
    assert list(bad)[0] in str(excinfo.value)


def test_tdz_names_a_compiled_operator(heisenberg4):
    """TDZ rebuilds its operators from their terms, so a toMPO() one
    must be refused with the message that says so, not die inside
    toMPO on a missing to_terms"""
    sc = heisenberg4
    B = sc.toMPO(sc.Sz[0])
    with pytest.raises(TypeError) as excinfo:
        sc.get_dynamical_correlator(mode="DMRG", submode="TDZ", name=[B, B],
                                    **BASE)
    assert "toMPO" in str(excinfo.value)
    assert "TDZ" in str(excinfo.value)


@pytest.mark.parametrize("submode", ["TD", "TDZ"])
@pytest.mark.parametrize("i,j", [(1, 1), (0, 2)])
def test_lower_level_route_honours_the_sites(heisenberg4, submode, i, j):
    """get_dynamical_correlator_MB(name="ZZ", i=, j=) is the explicit
    [Sz_i, Sz_j] pair; it used to be [Sz_0, Sz_0] bit for bit"""
    sc = heisenberg4
    y = np.asarray(sc.get_dynamical_correlator_MB(
        submode=submode, name="ZZ", i=i, j=j, **BASE)[1])
    y_ij = np.asarray(sc.get_dynamical_correlator(
        submode=submode, name=[sc.Sz[i], sc.Sz[j]], **BASE)[1])
    y_00 = np.asarray(sc.get_dynamical_correlator(
        submode=submode, name=[sc.Sz[0], sc.Sz[0]], **BASE)[1])
    assert np.max(np.abs(y - y_ij)) < 1e-12
    assert np.max(np.abs(y - y_00)) > 1e-3
