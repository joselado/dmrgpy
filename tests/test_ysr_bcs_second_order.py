"""Second-order spin-excitation dI/dV with superconducting electrodes,
ysrtk/bcsperturbative.py."""
import numpy as np
import pytest

from dmrgpy import spinchain, wilsonchain, ysr
from dmrgpy.kondospectrumtk import conductance
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
from dmrgpy.ysrtk import bcsperturbative as bp


def _spin1(T, kB=8.617333262e-5):
    sc = spinchain.Spin_Chain(["S=1"])
    sc.set_hamiltonian(1e-3*sc.Sz[0]*sc.Sz[0] + 0.3e-3*sc.Sx[0]
                       + 0.2e-3*sc.Sz[0])
    return KondoSpectrum(sc, 0, T, kB=kB)


@pytest.mark.parametrize("T", [0.0, 1.0, 4.0])
def test_no_gap_is_second_order_dIdV(T):
    ks = _spin1(T)
    eV = np.linspace(-4e-3, 4e-3, 401)
    ref = conductance.second_order_dIdV(ks, eV, U=0.3)
    got = bp.second_order_dIdV_bcs(ks, eV, U=0.3)
    assert np.max(np.abs(got - ref)) < 1e-13*np.max(ref)
    # the energy integral on the grid, without the closed-form normal part
    de = ks.kB*T/20 if T > 0 else 1e-6
    K = bp.bcs_kernel(ks, eV, de=de, force_grid=True, subtract_normal=False)
    grid = bp.second_order_dIdV_bcs(ks, eV, U=0.3, kernel=K)
    tol = 1e-10 if T == 0 else 5e-5
    assert np.max(np.abs(grid - ref)) < tol*np.max(ref)


def test_no_gap_current_is_ohmic_plus_steps():
    ks = _spin1(1.0)
    eV = np.linspace(-4e-3, 4e-3, 201)
    I = bp.second_order_current_bcs(ks, eV, U=0.3)
    G = conductance.second_order_dIdV(ks, eV, U=0.3)
    # trapezoid of G from 0, against I(V) - I(0)
    Iint = np.concatenate([[0], np.cumsum(0.5*(G[1:] + G[:-1])*np.diff(eV))])
    Iint = Iint - np.interp(0.0, eV, Iint)
    assert np.max(np.abs(I - Iint)) < 1e-3*np.max(np.abs(I))


def test_normal_tip_on_superconductor_is_the_dynes_dos():
    Delta, g = 1.0, 0.01
    K = bp.BCSKernel(3.5, 0.0, rho_s=lambda w: bp.dynes_dos(w, Delta, g),
                     de=g/20)
    x = np.linspace(-3, 3, 601) + 1e-4
    ok = np.abs(x) > 1e-3
    ref = bp.dynes_dos(x, Delta, g)*(x > 0)
    assert np.max(np.abs(K.conductance_ts(x) - ref)[ok]) < 1e-3*np.max(ref)
    ref = -bp.dynes_dos(x, Delta, g)*(x < 0)
    assert np.max(np.abs(K.conductance_st(x) - ref)[ok]) < 1e-3*np.max(-ref)


@pytest.mark.parametrize("Dt", [1.0, 0.6, 0.3])
def test_sis_current_jump(Dt):
    """At T=0 two BCS gaps give a jump pi sqrt(D1 D2)/2 of the current at
    D1 + D2 (in units of the normal-state current); the Dynes broadening
    makes it a Lorentzian of half-width g1 + g2."""
    g = 5e-4
    K = bp.BCSKernel(2.5, 0.0, rho_t=lambda w: bp.dynes_dos(w, Dt, g),
                     rho_s=lambda w: bp.dynes_dos(w, 1.0, g), de=g/20)
    x = np.linspace(1 + Dt - 0.05, 1 + Dt + 0.05, 40001)
    G = K.conductance_ts(x)
    weight = np.max(G)*np.pi*2*g
    assert weight == pytest.approx(np.pi/2*np.sqrt(Dt), rel=3e-3)
    assert abs(x[np.argmax(G)] - 1 - Dt) < 1e-6


def test_superconducting_tip_matches_tip_convolution():
    Delta, Dt, g = 1.0, 0.6, 0.01
    es = np.linspace(-6, 6, 24001)
    A = bp.dynes_dos(es, Delta, g)
    Vs = np.linspace(-3, 3, 121) + 1e-4
    _, ref = ysr.tip_convolution(es, A, Dt, gamma=g, Vs=Vs)
    K = bp.BCSKernel(3.5, 0.0, rho_s=lambda w: bp.dynes_dos(w, Delta, g),
                     rho_t=lambda w: bp.dynes_dos(w, Dt, g), de=g/20)
    got = K.conductance_ts(Vs) - K.conductance_st(Vs)
    assert np.max(np.abs(got - ref)) < 5e-4*np.max(ref)


def _brute_force_current(ks, eV, Ds, Dt, g, dxi, U, W=5.0):
    """Golden rule over the explicit Bogoliubov quasiparticles of a
    discretized BCS band in tip and sample: removing an electron destroys
    a quasiparticle (u^2 f(E), energy +E) or creates one (v^2 (1-f(E)),
    energy -E), adding one creates (u^2 (1-f(E)), +E) or destroys
    (v^2 f(E), -E). Every pair of events is one delta function of energy
    conservation, broadened into a Lorentzian of half-width 2g."""
    kT = ks.kB*ks.T
    def fermi(e):
        return 1.0/(np.exp(e/kT) + 1.0)
    def events(Delta, remove):
        xi = np.arange(-W + dxi/2, W, dxi)
        E = np.sqrt(xi**2 + Delta**2)
        u2, v2 = 0.5*(1 + xi/E), 0.5*(1 - xi/E)
        occ = fermi(E) if remove else 1.0 - fermi(E)
        return (np.concatenate([E, -E]),
                dxi*np.concatenate([u2*occ, v2*(1.0 - occ)]))
    b = g/5
    edges = np.arange(-2*W - 2, 2*W + 2 + b, b)
    cent = 0.5*(edges[1:] + edges[:-1])
    def pairs(e1, w1, e2, w2): # distribution of e2 - e1
        h = np.zeros(len(cent))
        for s in range(0, len(e1), 400):
            D = (e2[None, :] - e1[s:s+400, None]).ravel()
            Wt = (w1[s:s+400, None]*w2[None, :]).ravel()
            idx = np.floor((D - edges[0])/b).astype(int)
            h += np.bincount(idx, weights=Wt, minlength=len(cent))
        return h
    # t->s: eV = w_sample - w_tip + eps; s->t: eV = w_sample - w_tip - eps
    Hts = pairs(*events(Dt, True), *events(Ds, False))
    Hst = pairs(*events(Dt, False), *events(Ds, True))
    G2 = 2*g
    def rate(H, shift):
        d = eV[:, None] - cent[None, :] - shift
        return np.sum(H[None, :]*(G2/np.pi)/(d**2 + G2**2), axis=1)
    w, eps = bp._weights(ks)
    I = 4*U**2*(rate(Hts, 0.0) - rate(Hst, 0.0))
    for i in range(w.shape[0]):
        for f in range(w.shape[1]):
            if w[i, f] > 0:
                I = I + w[i, f]*(rate(Hts, eps[i, f]) - rate(Hst, -eps[i, f]))
    return 2*np.pi*I


def test_current_against_bogoliubov_quasiparticles():
    """Unequal gaps at kT = 0.15 Delta_s, so thermal quasiparticles take
    part. The two models differ at O(gamma): the explicit sum puts the
    occupation at the level energy, the Dynes form at the broadened one.
    So the difference must halve with gamma."""
    sc = spinchain.Spin_Chain(["S=1/2"])
    sc.set_hamiltonian(0.3*sc.Sz[0] + 0.1*sc.Sx[0])
    ks = KondoSpectrum(sc, 0, 0.15, kB=1.0)
    eV = np.linspace(-2.5, 2.5, 126)
    diffs = []
    for g, dxi in ((0.02, 0.006), (0.01, 0.003)):
        I = bp.second_order_current_bcs(ks, eV, Delta_s=1.0, Delta_t=0.6,
                                        gamma=g, U=0.2)
        Ib = _brute_force_current(ks, eV, 1.0, 0.6, g, dxi, U=0.2)
        diffs.append(np.max(np.abs(I - Ib))/np.max(np.abs(I)))
    assert diffs[1] < 5e-3
    assert diffs[1] < 0.6*diffs[0]


def _fe_oep_cl():
    """S=5/2 with D = 0.7 meV on Pb with a Pb tip, arXiv:1405.6984."""
    sc = spinchain.Spin_Chain(["S=5/2"])
    sc.set_hamiltonian(0.7e-3*sc.Sz[0]*sc.Sz[0])
    return sc


def _peaks(x, y, frac=1e-3):
    return np.array([x[i] for i in range(1, len(x) - 1)
                     if y[i] > y[i-1] and y[i] > y[i+1] and y[i] > frac*y.max()])


def test_fe_oep_cl_peaks():
    Delta = 1.35e-3
    eV = np.linspace(0.5e-3, 6e-3, 1101)
    _, G = bp.get_spin_excitation_spectrum(_fe_oep_cl(), eV, T=1.2,
            Delta_s=Delta, Delta_t=Delta, gamma=1e-2*Delta, U=1.15)
    p = _peaks(eV, G, 1e-2)
    assert np.allclose(p, [2*Delta, 2*Delta + 1.4e-3], atol=6e-6)
    ratio = np.interp(2*Delta + 1.4e-3, eV, G)/np.interp(2*Delta, eV, G)
    assert ratio == pytest.approx(0.441, abs=0.01)
    # at 4.8 K the state 3/2 is populated: 3/2 -> 5/2 at 2 Delta + 2.8 meV,
    # and equal gaps give the thermal peak at zero bias
    eV = np.linspace(-0.2e-3, 6e-3, 1241)
    _, G = bp.get_spin_excitation_spectrum(_fe_oep_cl(), eV, T=4.8,
            Delta_s=Delta, Delta_t=Delta, gamma=1e-2*Delta, U=1.15)
    p = _peaks(eV, G, 1e-2)
    for e in (0.0, 2*Delta, 2*Delta + 1.4e-3, 2*Delta + 2.8e-3):
        assert np.min(np.abs(p - e)) < 6e-6


def test_gap_needs_gamma():
    ks = _spin1(0.0)
    with pytest.raises(ValueError):
        bp.second_order_dIdV_bcs(ks, np.linspace(-1e-3, 1e-3, 11),
                                 Delta_s=1e-3)


def _spin1_rates(kT):
    sc = spinchain.Spin_Chain(["S=1"])
    sc.set_hamiltonian(0.5*sc.Sz[0]*sc.Sz[0])
    ks = KondoSpectrum(sc, 0, kT, kB=1.0)
    M2 = conductance._spin_matrix_elements_squared(ks)
    eps = ks.e[:, None] - ks.e[None, :]
    return ks, M2, eps, (M2 > 0.1) & (eps < -0.1)


def test_substrate_relaxation_rates():
    ks, M2, eps, down = _spin1_rates(0.0)
    R = bp.substrate_relaxation_rates(ks, 0.1)
    ref = np.pi*0.01*M2*np.maximum(-eps, 0.0)
    np.fill_diagonal(ref, 0.0)
    assert np.max(np.abs(R - ref)) < 1e-14
    # eps = 0.5 below 2 Delta = 0.6: only the Dynes in-gap states relax
    # it, so the rate goes to zero linearly in gamma
    r = [np.max(bp.substrate_relaxation_rates(ks, 0.1, Delta_s=0.3,
                gamma_s=g)[down]/ref[down]) for g in (1e-3, 1e-4)]
    assert r[0] < 1e-2
    assert r[1] == pytest.approx(r[0]/10, rel=2e-2)


def _bcs_rate_ratio(Delta, kT, x):
    """BCS rate of a spin-flip vertex absorbing x > 0 over the normal one,
    int f(w)(1 - f(w + x)) rho rho' (1 + Delta^2/(w w')) with signed
    energies, by direct quadrature: quasiparticle scattering above the gap
    (twice, by particle-hole symmetry) and pair creation."""
    from scipy.integrate import quad
    if kT > 0:
        f = lambda e: 0.5*(1 - np.tanh(e/(2*kT)))
    else:
        f = lambda e: float(e < 0)
    rho = lambda e: abs(e)/np.sqrt(e*e - Delta*Delta)
    tot = 0.0
    if kT > 0: # w = Delta cosh t
        def g(t):
            w = Delta*np.cosh(t)
            return (Delta*np.cosh(t)*rho(w + x)*(1 + Delta**2/(w*(w + x)))
                    *f(w)*(1 - f(w + x)))
        tot += 2*quad(g, 0, np.arccosh(1 + 60*kT/Delta) + 1)[0]
    if x > 2*Delta: # pair creation, E = -w in (Delta, x - Delta)
        g = lambda E: ((E*(x - E) - Delta**2)
                       /np.sqrt((E + Delta)*(x - E + Delta))
                       *f(-E)*(1 - f(x - E)))
        tot += quad(g, Delta, x - Delta, weight="alg", wvar=(-0.5, -0.5))[0]
    return tot/(x/(-np.expm1(-x/kT)) if kT > 0 else x)


@pytest.mark.parametrize("kT,Delta", [(0.0, 0.2), (0.0, 0.1), (0.1, 0.3),
                                      (0.1, 0.2)])
def test_relaxation_coherence_factor(kT, Delta):
    """The spin density is odd under time reversal: pair creation carries
    (1 - Delta^2/(E E')), so the rate at T=0 starts continuously at
    2 Delta, and quasiparticle scattering (1 + Delta^2/(E E')), the
    Hebel-Slichter factor. The Dynes rate converges to the BCS
    quadrature linearly in gamma."""
    ks, M2, eps, down = _spin1_rates(kT)
    R0 = bp.substrate_relaxation_rates(ks, 0.1, de=1e-4 if kT == 0 else None)
    ref = _bcs_rate_ratio(Delta, kT, 0.5)
    d = [np.max(bp.substrate_relaxation_rates(ks, 0.1, Delta_s=Delta,
                gamma_s=g)[down]/R0[down]) - ref for g in (4e-4, 2e-4)]
    assert abs(d[1]) < 2e-3 # in units of the normal-metal rate
    assert d[1] == pytest.approx(d[0]/2, rel=0.1)


def test_kernel_reuse_is_checked():
    sc = spinchain.Spin_Chain(["S=1/2"])
    sc.set_hamiltonian(0.2e-3*sc.Sz[0])
    Delta, g = 1.35e-3, 1.35e-5
    ks = KondoSpectrum(sc, 0, 1.2)
    K = bp.bcs_kernel(ks, np.linspace(-2e-3, 2e-3, 101), Delta_s=Delta,
                      gamma_s=g)
    eV = np.linspace(-2e-3, 2e-3, 51)
    assert np.allclose(bp.second_order_dIdV_bcs(ks, eV, kernel=K),
            bp.second_order_dIdV_bcs(ks, eV, Delta_s=Delta, gamma=g))
    with pytest.raises(ValueError, match="covers"):
        bp.second_order_dIdV_bcs(ks, np.linspace(1e-3, 4e-3, 31), kernel=K)
    with pytest.raises(ValueError, match="kT"):
        bp.second_order_dIdV_bcs(KondoSpectrum(sc, 0, 4.8), eV, kernel=K)


def test_callable_dos_without_any_scale_needs_de():
    ks = _spin1(0.0)
    flat = lambda w: np.ones_like(w)
    with pytest.raises(ValueError, match="de"):
        bp.second_order_dIdV_bcs(ks, np.linspace(-1e-3, 1e-3, 11),
                                 rho_s=flat)
    with pytest.raises(ValueError, match="de"):
        bp.substrate_relaxation_rates(ks, 0.1, rho_s=flat)
    with pytest.raises(ValueError, match="pair_s"):
        bp.substrate_relaxation_rates(ks, 0.1, Delta_s=1e-3, rho_s=flat,
                                      de=1e-6)


def test_many_body_spectrum_at_zero_exchange():
    """At J=0 the composite-operator spectrum of ysr.get_stm_spectrum is
    the golden rule with the local density of states of the first site of
    the BCS Wilson chain: the same curve, up to the Lorentzian tails the
    golden rule cuts at the Fermi level, and the factor 4 pi of the two
    normalizations (U=0)."""
    n, Delta, B, Lam, delta = 3, 0.1, 0.15, 4.0, 0.005
    es = np.linspace(-0.9, 0.9, 361)
    e, t = wilsonchain.wilson_chain(n, Lambda=Lam, band=(-1, 1))
    h = np.diag(e) + np.diag(t, 1) + np.diag(t, -1)
    E, V = np.linalg.eigh(np.block([[h, Delta*np.eye(n)],
                                    [Delta*np.eye(n), -h]]))
    w0 = np.abs(V[0, :])**2
    def rho(w):
        return np.sum(w0[None, :]*(delta/np.pi)
                      /((np.asarray(w)[:, None] - E[None, :])**2 + delta**2),
                      axis=1)
    sc = spinchain.Spin_Chain(["S=1/2"])
    sc.set_hamiltonian(B*sc.Sz[0])
    ks = KondoSpectrum(sc, 0, 0.0)
    Gp = bp.second_order_dIdV_bcs(ks, es, rho_s=rho, de=delta/20)/(4*np.pi)
    fc = ysr.kondo_chain(n, J=0.0, Delta=Delta, Lambda=Lam, B=B)
    _, Gm = ysr.get_stm_spectrum(fc, es, delta=delta, mode="ED",
                                 submode="CVM")
    assert np.max(np.abs(Gm - Gp)) < 3e-3*np.max(Gm)
    assert np.allclose(_peaks(es, Gm, 0.02), _peaks(es, Gp, 0.02))
