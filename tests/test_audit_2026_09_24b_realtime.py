"""Regressions for the real-time cluster of the 2026-09-24b hole hunt
(docs/audit_2026_09_24b_hole_hunt.md, findings 7 to 10).

Every test here is held to an anchor built outside dmrgpy (a closed form,
or a matrix exponential written with the sign of the Schrodinger equation
fixed by hand), never to another backend, since findings 9 and 10 moved
ED and DMRG at the same time and agreement between them was what hid
both.

- Finding 9. `edtk/timedependent.evolution_ABC`, behind
  `evolve_and_measure(mode="ED")` and `evolution_ABA(mode="ED")`, evolved
  with e^{+iHt} and so returned the time-reversed trajectory. Pinned by
  Larmor precession, <Sy_0>(t) = +sin(Bt)/2 and <S+_0>(t) = e^{+iBt}/2,
  and by one fermion on a ring with a flux, which has to go round the
  way the flux says.
- Finding 10. DMRG `evolve_and_measure` conjugated its result, so it
  returned <psi(t)|O^dagger|psi(t)>. Pinned by an eigenstate, where the
  answer is vev(O) at every t with no propagator involved, on every
  backend and integrator, ED included.
- Finding 8. `sxt_to_skomega` returned every infinite-chain S(k,omega)
  mirrored in omega. Pinned on a chiral free-fermion ring, the one kind
  of model where conjugating before the spatial sum (the wrong fix) and
  after it (the right one) differ, by the momentum label.
- Finding 7. `sxt_to_skomega` was hard-wired to the FFT-plus-
  interpolation stage. Pinned against the closed form of the damped
  trapezoid sum at every window, the public defaults included.
"""
import numpy as np
import pytest
import scipy.linalg as sla

from dmrgpy import cppext, fermionchain, spinchain, timedependent

from _helpers import julia_available


def _row(mode, version, *rest):
    """pytest.param(mode, version, *rest), skipped when a compiled backend
    is not there. The julia_live check is deferred to the test body, so
    that collecting this file (and `-k "not julia_live"`) never boots a
    Julia session."""
    ident = "-".join([mode, "v%d" % version if version in (2, 3)
                      else str(version)] + [str(r) for r in rest])
    marks = ()
    if version in (2, 3):
        marks = pytest.mark.skipif(not cppext.available(version),
            reason="needs the compiled itensor_version=%d extension"
                   % version)
    return pytest.param(mode, version, *rest, id=ident, marks=marks)


BACKENDS = [_row("ED", "python"), _row("DMRG", "python"), _row("DMRG", 3)]


def _require_julia(version):
    if version == "julia_live":
        ok, reason = julia_available()
        if not ok:
            pytest.skip("requires a working juliacall/Julia toolchain: %s"
                        % reason)


def _spin_chain(version, n, tevol_method="TDVP"):
    sc = spinchain.Spin_Chain([2] * n, itensor_version=version)
    sc.maxm, sc.nsweeps = 10, 10
    sc.tevol_method = tevol_method
    return sc


# ------------------------------------------------------------- finding 9

B_FIELD, NT9, DT9 = 1.0, 40, 0.1
TS9 = DT9 * np.arange(NT9)


@pytest.mark.parametrize("mode,version", BACKENDS)
def test_larmor_precession_runs_forward(mode, version):
    """H = B sum Sz from the +x product state: dSy/dt = i[H,Sy] = B Sx, so
    <Sx_0>(t) = cos(Bt)/2, <Sy_0>(t) = +sin(Bt)/2 and <S+_0>(t) =
    e^{+iBt}/2. The quench is the one the examples write: ground state of
    -sum Sx, then set_hamiltonian(B sum Sz). ED used to return
    <Sy_0>(t=1) = -0.4207 against +0.4207, and <S+_0> agreed with DMRG
    only because the two were wrong in the same way (DMRG by finding 10's
    conjugation)."""
    n = 3
    sc = _spin_chain(version, n)
    sc.set_hamiltonian(sum(-sc.Sx[i] for i in range(n)))
    wf = sc.get_gs(mode=mode)
    sc.set_hamiltonian(sum(B_FIELD * sc.Sz[i] for i in range(n)))
    cases = ((sc.Sx[0], 0.5 * np.cos(B_FIELD * TS9)),
             (sc.Sy[0], 0.5 * np.sin(B_FIELD * TS9)),
             (sc.Sx[0] + 1j * sc.Sy[0], 0.5 * np.exp(1j * B_FIELD * TS9)))
    for op, exact in cases:
        _t, y = timedependent.evolve_and_measure(sc, operator=op, nt=NT9,
                                                  dt=DT9, wf=wf, mode=mode)
        # solve_ivp's own tolerance is ~1e-8 here, DMRG ~1e-14; the
        # backward trajectory is off by up to 1.0
        assert np.max(np.abs(np.asarray(y) - exact)) < 1e-5


def _flux_ring(version, n=3, hop=1.0, phi=np.pi / 6, mu=3.0):
    fc = fermionchain.Fermionic_Chain(n, itensor_version=version)
    h = 0
    for j in range(n):
        h = h - hop * np.exp(1j * phi) * fc.Cdag[(j + 1) % n] * fc.C[j]
        h = h - hop * np.exp(-1j * phi) * fc.Cdag[j] * fc.C[(j + 1) % n]
    for j in range(n):
        h = h + mu * fc.N[j]
    fc.set_hamiltonian(h)
    fc.maxm, fc.nsweeps = 20, 10
    fc.tevol_method = "TDVP"
    hsp = np.diag([mu] * n).astype(complex)
    for j in range(n):
        hsp[(j + 1) % n, j] += -hop * np.exp(1j * phi)
        hsp[j, (j + 1) % n] += -hop * np.exp(-1j * phi)
    return fc, hsp


@pytest.mark.parametrize("mode,version", BACKENDS)
def test_flux_ring_circulates_forward(mode, version):
    """One fermion on a 3-site ring with flux pi/6, mu = 3 > 2t so the
    vacuum is the unique ground state and Cdag_0|vac> is one particle on
    site 0. In the one-particle sector the Jordan-Wigner strings only
    cross empty sites, so <N_k>(t) = |(e^{-i h_sp t} e_0)_k|^2 exactly.
    The flux breaks time reversal and time reversal swaps N_1 and N_2,
    so a backward run is off by up to 1.0 on a real observable, where
    conjugating the result could not help: ED used to return
    <N_1>(t=1.5) = 0.1024 against 0.8413."""
    fc, hsp = _flux_ring(version)
    e0 = np.zeros(3, complex)
    e0[0] = 1.0
    fw = np.array([np.abs(sla.expm(-1j * hsp * t) @ e0) ** 2 for t in TS9])
    for k in (1, 2):
        _t, y = timedependent.evolution_ABA(fc, A=fc.Cdag[0], B=fc.N[k],
                                             mode=mode, nt=NT9, dt=DT9)
        assert np.max(np.abs(np.asarray(y) - fw[:, k])) < 1e-4
    assert np.max(np.abs(fw[:, 1] - fw[:, 2])) > 0.9  # the anchor sees it


# ------------------------------------------------------------ finding 10

NT10, DT10 = 6, 0.1
EIGEN_ROWS = (
    [_row("ED", "python", "TDVP")]
    + [_row("DMRG", "python", m)
       for m in ("TDVP", "TDVP_GSE", "TEBD", "AUTO", "MPO")]
    + [_row("DMRG", 3, m) for m in ("TDVP", "TDVP_GSE", "TEBD", "AUTO", "MPO")]
    + [_row("DMRG", 2, "MPO")]
    + [_row("DMRG", "julia_live", m) for m in ("TDVP", "TEBD")]
)


def _plus_x_chain(mode, version, tevol_method, n=3):
    _require_julia(version)
    sc = _spin_chain(version, n, tevol_method)
    sc.set_hamiltonian(sum(-sc.Sx[i] for i in range(n)))
    wf = sc.get_gs(mode=mode)
    return sc, wf


@pytest.mark.parametrize("mode,version,tevol_method", EIGEN_ROWS)
def test_evolve_and_measure_on_an_eigenstate_is_vev(mode, version,
                                                     tevol_method):
    """H = -sum Sx has the +x product state as its unique ground state, an
    eigenstate, so <O>(t) is constant and equal to vev(O) whatever the
    propagator. For O = Sz_0 + i*Sx_0 that is exactly +0.5i. DMRG used to
    return -0.5i at every t while vev() returned +0.5i on the same state.
    The t=0 value is measured before any step, so it is vev(O) to
    rounding on every integrator; the trajectory is held to 1e-3, which
    lets through the MPO-Taylor stepper's own 6e-4."""
    sc, wf = _plus_x_chain(mode, version, tevol_method)
    op = sc.Sz[0] + 1j * sc.Sx[0]
    ref = complex(sc.vev(op, mode=mode))
    assert ref == pytest.approx(0.5j, abs=1e-8)
    _t, y = timedependent.evolve_and_measure(sc, operator=op, nt=NT10,
                                              dt=DT10, wf=wf, mode=mode)
    y = np.asarray(y)
    assert y[0] == pytest.approx(ref, abs=1e-8)
    assert np.max(np.abs(y - 0.5j)) < 1e-3


@pytest.mark.parametrize("mode,version", BACKENDS[1:])
def test_evolve_and_measure_return_wf_branch_is_not_conjugated(mode, version):
    """The return_wf=True branch had its own copy of the conjugating
    return."""
    sc, wf = _plus_x_chain(mode, version, "TDVP")
    op = sc.Sz[0] + 1j * sc.Sx[0]
    _t, y, wf_final = timedependent.evolve_and_measure(
        sc, operator=op, nt=NT10, dt=DT10, wf=wf, return_wf=True)
    assert np.max(np.abs(np.asarray(y) - 0.5j)) < 1e-8
    assert wf_final is not None


@pytest.mark.parametrize("mode,version", BACKENDS)
def test_evolution_aba_on_an_eigenstate_is_vev(mode, version):
    """evolution_ABA, the one consumer of evolve_and_measure_dmrg inside
    src/: A = 2*Sx_1 leaves the +x state unchanged, so <B>(t) with B =
    Sz_0 + i*Sx_0 is +0.5i at every t."""
    sc, wf = _plus_x_chain(mode, version, "TDVP")
    kw = dict(wf=wf) if mode == "DMRG" else {}
    _t, y = timedependent.evolution_ABA(sc, A=2 * sc.Sx[1],
                                         B=sc.Sz[0] + 1j * sc.Sx[0],
                                         nt=NT10, dt=DT10, mode=mode, **kw)
    assert np.max(np.abs(np.asarray(y) - 0.5j)) < 1e-6


# ------------------------------------------------------------- finding 8

def test_sxt_to_skomega_momentum_label_and_sign_on_a_chiral_ring():
    """A free-fermion ring with complex first-neighbour hopping t e^{i phi}
    (band minimum off k=0, Fermi sea off centre) and the pair (Cdag_x,
    C_0), no DMRG. S(x,t) is built exactly as the IBC window produces it,
    <0|A_x e^{-i(H-E0)t} B_0|0> = sum_q M_q(x) e^{-i D_q t}, with
    |q> = c_q|0>, D_q = -eps_q > 0 and M_q(x) = conj(v[x,q]) v[0,q] on
    the occupied levels. The reference is the house density of the
    momentum series, sum_x e^{-ikx} sum_q M_q(x) L_delta(w - D_q), lines
    at +D_q. Before the fix every line came out at -D_q (error 1.0 of the
    scale); conjugating each x series before the spatial sum instead
    puts every k's band at -k (error 1.0 again), and the two k with no
    occupied state, whose reference is zero, are what catch that."""
    N, hop, phi, mu = 40, 1.0, 0.7, -0.4
    h = np.zeros((N, N), dtype=complex)
    for x in range(N):
        h[(x + 1) % N, x] += -hop * np.exp(1j * phi)
        h[x, (x + 1) % N] += -hop * np.exp(-1j * phi)
    h -= mu * np.eye(N)
    eps, v = np.linalg.eigh(h)
    occ = eps < 0
    D = -eps[occ]
    xs = np.arange(N)
    M = np.conj(v[:, occ]) * v[0, occ][None, :]            # M[x, q]
    dt, nt, delta = 0.05, 1600, 0.2                           # delta*T = 16
    ts = dt * np.arange(nt)
    S = (M[None, :, :] * np.exp(-1j * D[None, None, :]
                                * ts[:, None, None])).sum(-1)
    es = np.linspace(-3.5, 3.5, 701)
    ms = np.array([-7, -5, -3, 3, 5, 7])
    ks = 2 * np.pi * ms / N
    _k, _e, Skw = timedependent.sxt_to_skomega(ts, xs, S, dt, ks=ks, es=es,
                                                delta=delta)
    lor = delta / np.pi / ((es[:, None] - D[None, :]) ** 2 + delta ** 2)
    ref = np.array([(lor * (np.exp(-1j * k * xs) @ M)[None, :]).sum(1)
                    for k in ks])
    scale = np.max(np.abs(ref))
    assert np.max(np.abs(Skw.real - ref.real)) / scale < 1e-4
    # the anchor is chiral: the reference at k and at -k differ by O(1)
    assert np.max(np.abs(ref[:3][::-1] - ref[3:])) / scale > 0.5


# ------------------------------------------------------------- finding 7

def _single_magnon(dt, nt, L=16):
    """S(x,t) = (1/L) sum_q e^{iqx} e^{-i eps(q) t}, eps(q) = 2 - 2cos q,
    the IBC-window sign, and its momentum grid."""
    xs = np.arange(L) - L // 2
    qs = 2 * np.pi * np.arange(L) / L
    qs = np.where(qs > np.pi, qs - 2 * np.pi, qs)
    ts = dt * np.arange(nt)
    S = np.exp(-1j * np.outer(ts, 2 - 2 * np.cos(qs))) @ (
        np.exp(1j * np.outer(qs, xs)) / L)
    ks = np.sort(qs)
    return ts, xs, S, ks, 2 - 2 * np.cos(ks)


@pytest.mark.parametrize("dt,nt,delta", [
    (0.1, 200, 0.05),     # td_dynamical_correlator's defaults, delta*T = 1
    (0.1, 200, 0.3),      # delta*T = 6 at the default dt and nt
    (0.1, 1200, 0.05),    # delta*T = 6 at the default delta
    (0.05, 40, 0.15),     # the example's heatmap call, delta*T = 0.29
])
def test_sxt_to_skomega_is_the_direct_sum_at_every_window(dt, nt, delta):
    """On the q grid the conjugated momentum series of a single magnon is
    exactly e^{+i eps(k) t}, whose damped trapezoid sum is a geometric
    series, (dt/pi)[sum_{j<n} r^j - (1 + r^{n-1})/2] with
    r = e^{(i(eps - w) - delta) dt}. The reduction used to interpolate an
    FFT of spacing 2*pi/(nt*dt) onto es whatever the window, which is
    2.1e-01 of the peak off this at delta*T = 6 and more at the short
    windows. Run on the default es (window [-1,10], 800 points), so the
    public defaults are covered too."""
    ts, xs, S, ks, eps_k = _single_magnon(dt, nt)
    _ks, es, Skw = timedependent.sxt_to_skomega(ts, xs, S, dt, ks=ks,
                                                 delta=delta)
    r = np.exp((1j * (eps_k[:, None] - es[None, :]) - delta) * dt)
    ref = dt / np.pi * ((1 - r ** nt) / (1 - r) - 0.5 * (1 + r ** (nt - 1)))
    assert np.max(np.abs(Skw - ref)) / np.max(np.abs(ref)) < 1e-11
