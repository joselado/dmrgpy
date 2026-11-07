"""Current through a YSR state with a superconducting tip
(ysrtk/ysrcurrent.py, Ruby et al., arXiv:1502.05048): the quadrature, the
linear regime against ysr.tip_convolution, the threshold values of the
paper's supplement, the peak ratios, the thermal peaks, the normal-tip
Andreev peak, and the rate parameters from a many-body chain."""
import warnings

import numpy as np
import pytest
from scipy.integrate import quad

from dmrgpy import wilsonchain, ysr
from dmrgpy.ysrtk import ysrcurrent as yc


def test_graded_mesh_agrees_with_quad():
    p = dict(eps=0.3, u2=0.2, v2=1.0, Delta_t=1.0, gamma_tip=1e-3,
             Gamma1=1e-3, Gamma2=2e-4, T=0.05, dynes=1e-3)
    for V in (1.3, -1.3, 0.7):
        Is, Ia = yc.current([V], **p)
        W = 10*(1 + 0.3 + abs(V))
        pts = sorted(np.clip([0.3, V - 1, V + 1, -V - 1, -V + 1], -W, W))
        def f(w, k):
            return yc._integrands(np.array([w]), V, 0.3, 0.2, 1.0, 1.0, 1e-3,
                                  1e-3, 2e-4, 0.05, 1e-3)[k][0]
        for k, I in ((0, Is[0]), (1, Ia[0])):
            q = 0.5*quad(f, -W, W, args=(k,), points=pts, limit=2000,
                         epsabs=1e-16, epsrel=1e-12)[0]
            assert I == pytest.approx(q, rel=1e-10)


def test_linear_regime_is_the_tip_convolution():
    """gamma_e << Gamma_1, Gamma_2 = 0, T = 0: Gamma_1/[(w-eps)^2 +
    Gamma_1^2/4] is 2 pi times a Lorentzian of half-width Gamma_1/2, so
    dI_s/dV = pi g d/dV of tip_convolution of u2 L(w-eps) + v2 L(w+eps)."""
    Dt, eps, u2, v2, G1, dy, g = 1.0, 0.3, 0.2, 1.0, 4e-3, 1e-2, 1e-8
    es = np.linspace(-4, 4, 160001)
    L = lambda x: (G1/2)/np.pi/(x*x + (G1/2)**2)
    A = u2*L(es - eps) + v2*L(es + eps)
    Vs = np.concatenate([np.linspace(1.2, 1.45, 21), -np.linspace(1.2, 1.45, 21)])
    ds, da = yc.didv(Vs, eps, u2, v2, Dt, g, G1, 0.0, dynes=dy, window=4.0)
    _, dc = ysr.tip_convolution(es, A, Dt, dy, Vs=Vs)
    dc = np.pi*g*dc
    assert np.max(np.abs(ds - dc)) < 2e-4*np.max(np.abs(ds))
    assert np.max(np.abs(da)) < 1e-5*np.max(np.abs(ds))


def test_threshold_values_of_the_supplement():
    """S68 (Gamma_1 >> omega_e), S71 and S70 (omega_e >> Gamma_1), at
    eV = Delta + eps, with v2 = 0 so that there is no Andreev current."""
    Dt, eps = 1.0, 0.3
    thr = Dt + eps
    ge, G1 = 1e-9, 1e-3 # regime (i)
    ds, _ = yc.didv([thr], eps, 1.0, 0.0, Dt, ge, G1, 0.0)
    assert ds[0] == pytest.approx(np.pi/np.sqrt(2)*ge*np.sqrt(Dt)/G1**1.5,
                                  rel=1e-4)
    ge, G1 = 1e-4, 1e-8 # regime (ii)
    we = 0.5*(ge*np.sqrt(Dt))**(2.0/3.0)
    Is, _ = yc.current([thr], eps, 1.0, 0.0, Dt, ge, G1, 0.0)
    ds, _ = yc.didv([thr], eps, 1.0, 0.0, Dt, ge, G1, 0.0, h=0.02*we)
    assert Is[0] == pytest.approx(np.pi*G1/3, rel=2e-3) # e Gamma_1/(3 hbar)
    assert ds[0] == pytest.approx(2*np.pi*G1/(9*we), rel=1e-3)


def test_peak_ratios_in_the_linear_regime():
    """alpha_+/alpha_- = u2/v2, alpha_+ beta_+ = alpha_- beta_-,
    beta_-/alpha_+ = Gamma_2/Gamma_1 (S57, S76, S80, S83)."""
    eps, u2, v2, G1, G2 = 0.3, 0.13, 1.0, 1e-3, 2.5e-4
    p = yc.peak_heights(eps, u2, v2, 1.0, 1e-8, G1, G2)
    h = {k: p[k]["total"] for k in p}
    assert h["alpha+"]/h["alpha-"] == pytest.approx(u2/v2, rel=5e-3)
    assert h["alpha+"]*h["beta+"]/(h["alpha-"]*h["beta-"]) == pytest.approx(1, abs=1e-3)
    assert h["beta-"]/h["alpha+"] == pytest.approx(G2/G1, rel=1e-3)
    assert p["alpha+"]["a"] < 1e-4*p["alpha+"]["s"]


def test_asymmetry_inverts_at_strong_tunneling():
    """The Andreev current takes over first at the peak with the smaller
    weight, and alpha_+/alpha_- inverts (Fig. 3 of the paper)."""
    eps, u2, v2, G1, G2 = 0.3, 0.13, 1.0, 1e-3, 2.5e-4
    lo = yc.peak_heights(eps, u2, v2, 1.0, 1e-8, G1, G2, peaks=("alpha+", "alpha-"))
    lo2 = yc.peak_heights(eps, u2, v2, 1.0, 1e-7, G1, G2, peaks=("alpha+",))
    hi = yc.peak_heights(eps, u2, v2, 1.0, 1e-2, G1, G2, peaks=("alpha+", "alpha-"))
    assert np.log10(lo2["alpha+"]["total"]/lo["alpha+"]["total"]) == pytest.approx(1, abs=1e-2)
    assert lo["alpha+"]["total"] < lo["alpha-"]["total"]
    assert hi["alpha+"]["total"] > 4*hi["alpha-"]["total"]
    assert hi["alpha+"]["a"] > 10*hi["alpha+"]["s"]


def test_thermal_peaks_follow_detailed_balance():
    """At T > 0 with Gamma_2 = Gamma_1 exp(-eps/T) the thermal peaks sit at
    +-(Delta_t - eps) with beta_-/alpha_+ = exp(-eps/T), up to the thermal
    quasiparticles of the tip, exp(-(Delta_t - eps)/T) of beta_-."""
    eps, T = 0.3, 0.1
    p = yc.peak_heights(eps, 0.13, 1.0, 1.0, 1e-7, 1e-3, None, T=T)
    assert p["beta-"]["total"]/p["alpha+"]["total"] == pytest.approx(np.exp(-eps/T), rel=2e-3)
    assert p["beta+"]["total"]/p["alpha-"]["total"] == pytest.approx(np.exp(-eps/T), rel=2e-3)
    assert p["beta+"]["V"] == pytest.approx(0.7, abs=1e-3)
    assert p["beta-"]["V"] == pytest.approx(-0.7, abs=1e-3)
    q = yc.peak_heights(eps, 0.13, 1.0, 1.0, 1e-7, 1e-3, None, T=0.0)
    assert q["beta-"]["total"] < 1e-6*q["alpha+"]["total"]


def test_normal_tip_resonant_andreev_peak():
    """Delta_t = 0, no relaxation: a Breit-Wigner Andreev peak at eV =
    +-eps of height 4 Gamma_e Gamma_h/(Gamma_e + Gamma_h)^2 in 2e^2/h, the
    same at both signs of the bias."""
    eps, u2, v2 = 0.3, 0.13, 1.0
    V = np.linspace(0.28, 0.32, 401)
    ds, da = yc.didv(V, eps, u2, v2, 0.0, 1e-2, 0.0, 0.0)
    _, dan = yc.didv(-V, eps, u2, v2, 0.0, 1e-2, 0.0, 0.0)
    assert np.max(np.abs(ds)) == 0.0
    assert np.max(da) == pytest.approx(4*u2*v2/(u2 + v2)**2, rel=1e-4)
    assert np.max(dan) == pytest.approx(np.max(da), rel=1e-10)


def test_gamma_tip_from_GN():
    for GN in (1e-6, 0.1, 1.0):
        g = yc.gamma_tip_from_GN(GN, nu_s=0.5)
        x = np.pi*0.5*g/2 # pi^2 nu_t nu_s t^2 with g = 2 pi t^2 nu_t
        assert 4*x/(1 + x) == pytest.approx(GN, rel=1e-12)
        g = yc.gamma_tip_from_GN(GN, nu_s=0.5, exact="landauer")
        x = np.pi*0.5*g/2
        assert x <= 1.0 + 1e-12
        assert 4*x/(1 + x)**2 == pytest.approx(GN, rel=1e-12)
    assert yc.gamma_tip_from_GN(1e-8) == pytest.approx(1e-8/np.pi, rel=1e-8)
    assert yc.gamma_tip_from_GN(1e-8, exact="landauer") == pytest.approx(1e-8/np.pi, rel=1e-8)
    assert yc.gamma_tip_from_GN(1e-2, exact=False) == pytest.approx(1e-2/np.pi, rel=1e-14)
    r = yc.gamma_tip_from_GN(0.1, exact="landauer")/yc.gamma_tip_from_GN(0.1)
    assert r == pytest.approx(1.027, abs=1e-3)
    with pytest.raises(ValueError):
        yc.gamma_tip_from_GN(4.0)
    with pytest.raises(ValueError):
        yc.gamma_tip_from_GN(1.5, exact="landauer")
    with pytest.raises(ValueError):
        yc.gamma_tip_from_GN(0.1, exact="paper")


def test_andreev_asymmetry_and_regime_iii():
    """With a small Gamma_1 the ratio alpha_+/alpha_- approaches the S57
    value, near (v/u)^(10/3) = 30, and then falls once the hole broadening
    sqrt(omega_h^3/eps) reaches omega_e (regime (iii), S65, S85); with
    Gamma_1 = 1e-3 there is no window and it stays near 5."""
    eps, u2, v2 = 0.3, 0.13, 1.0
    r = {}
    for G1, g in ((1e-7, 1e-5), (1e-7, 1e-3), (1e-3, 1e-2)):
        p = yc.peak_heights(eps, u2, v2, 1.0, g, G1, 0.0, peaks=("alpha+", "alpha-"))
        r[G1, g] = p["alpha+"]["total"]/p["alpha-"]["total"]
    assert r[1e-7, 1e-5] == pytest.approx(23.47, rel=1e-3)
    assert r[1e-7, 1e-3] == pytest.approx(12.75, rel=1e-3)
    assert r[1e-3, 1e-2] == pytest.approx(4.987, rel=1e-3)


def _split_poles(e1, e2, gap=0.1):
    """Electron poles at e1, e2 (weights 0.04, 0.02), hole poles at -e1,
    -e2 (0.03, 0.01), and a continuum at +-0.5 that makes the total weight
    2 (norm 1)."""
    om = np.array([e1, e2, -e1, -e2, 0.5, -0.5])
    wt = np.array([0.04, 0.02, 0.03, 0.01, 0.9, 1.0])
    return dict(omega=om, weight=wt)


def test_rate_parameters_sums_a_split_doublet():
    """A doublet split by a ground state of mixed parity (the subgap poles
    of an 8-site chain near the transition, at 0.01765 and 0.01833 Delta)
    is summed, with a warning, rather than halved."""
    with pytest.warns(RuntimeWarning, match="split"):
        r = yc.rate_parameters(_split_poles(1.764956e-3, 1.833468e-3), 0.1)
    assert r["u2"] == pytest.approx(0.06, rel=1e-12)
    assert r["v2"] == pytest.approx(0.04, rel=1e-12)
    assert r["norm"] == pytest.approx(1.0, rel=1e-12)
    assert len(r["others"]["omega"]) == 0
    # split beyond rtol*|gap|: only the strongest pole, and a warning on the rest
    with pytest.warns(RuntimeWarning, match="away from"):
        r = yc.rate_parameters(_split_poles(1.2e-3, 1.6e-3), 0.1)
    assert r["u2"] == pytest.approx(0.04, rel=1e-12)
    assert len(r["others"]["omega"]) == 2
    # a parity mixture flagged by the elastic weight of stm_poles
    p = _split_poles(1.9e-3, 1.9e-3)
    p["elastic"] = 1e-3
    with pytest.warns(RuntimeWarning, match="elastic"):
        yc.rate_parameters(p, 0.1)
    # one converged doublet: no warning, and the check against E_ysr
    p = _split_poles(1.9e-3, 1.9e-3 + 1e-12)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        r = yc.rate_parameters(p, 0.1, E_ysr=-1.9e-3)
    assert r["u2"] == pytest.approx(0.06, rel=1e-12)
    with pytest.raises(ValueError):
        yc.rate_parameters(p, 0.1, E_ysr=1.7e-3)


def test_negative_eps_raises():
    with pytest.raises(ValueError):
        yc.current([1.0], -0.1, 1.0, 1.0, 1.0, 1e-3, 1e-3)


def test_chain_weight_times_A_Lambda_is_the_continuum_residue():
    """Classical spin on the BdG Wilson chain: A_Lambda times the weight of
    the subgap pole at the first site is the residue u^2, v^2 of S32."""
    D, Lambda, n, Delta = 1.0, 2.0, 40, 1e-3
    A = ysr.lambda_factor(Lambda)
    nu0 = 1/(2*D)
    for JS, V in ((0.4, 0.0), (0.6, 0.2)):
        a, b = np.pi*nu0*JS, np.pi*nu0*V
        den = ((1 - a*a + b*b)**2 + 4*a*a)**1.5
        u2 = 2*a*np.pi*nu0*Delta*(1 + (a + b)**2)/den
        v2 = 2*a*np.pi*nu0*Delta*(1 + (a - b)**2)/den
        eps, t = wilsonchain.wilson_chain(n, Lambda=Lambda, band=(-D, D))
        h = np.diag(eps) + np.diag(t, 1) + np.diag(t, -1)
        hu, hd = h.copy(), h.copy()
        hu[0, 0] += A*(V - JS)
        hd[0, 0] += A*(V + JS)
        H = np.block([[hu, Delta*np.eye(n)], [Delta*np.eye(n), -hd]])
        E, W = np.linalg.eigh(H)
        k = np.flatnonzero((E > 0) & (E < Delta))[0]
        assert A*W[0, k]**2 == pytest.approx(u2, rel=1e-2)
        assert A*W[n, k]**2 == pytest.approx(v2, rel=1e-2)


def test_rate_parameters_from_a_kondo_chain():
    """eps, u2, v2 from the TASK poles of a 5-site chain: eps is
    |E_ysr|, the bare electron has u2 = v2 (particle-hole symmetry) and
    total weight 1 per spin, the composite vertex with U = 0.3 has
    u2/v2 = 1.557."""
    from dmrgpy.ysrtk import zaverage
    Lambda, Delta = 4.0, 0.1
    fc = ysr.kondo_chain(5, J=0.8*ysr.lambda_factor(Lambda), Delta=Delta,
                         Lambda=Lambda)
    fc.maxm = 60
    fc.nsweeps = 15
    r = ysr.subgap_energies(fc, mode="ED")
    s = fc.channels[0][0]
    pb = zaverage.stm_poles(fc, operators=[fc.Cdagup[s], fc.Cdagdn[s]], nkry=100)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        rb = yc.rate_parameters(pb, Delta, Lambda=Lambda, E_ysr=r["E_ysr"])
    assert rb["eps"] == pytest.approx(abs(r["E_ysr"]), abs=1e-8)
    assert rb["spread"] < 1e-6*Delta
    assert rb["norm"] == pytest.approx(1.0, abs=1e-10)
    assert rb["u2"] == pytest.approx(rb["v2"], rel=1e-6)
    assert rb["u2"] == pytest.approx(ysr.lambda_factor(Lambda)*rb["We"], rel=1e-12)
    px = zaverage.stm_poles(fc, U=0.3, nkry=100)
    rx = yc.rate_parameters(px, Delta, Lambda=Lambda)
    assert rx["eps"] == pytest.approx(rb["eps"], abs=1e-8)
    assert rx["u2"]/rx["v2"] == pytest.approx(1.557234, rel=1e-4)
    sub = dict(omega=px["omega"][np.abs(px["omega"]) < Delta],
               weight=px["weight"][np.abs(px["omega"]) < Delta])
    with pytest.raises(ValueError):
        yc.rate_parameters(sub, Delta)
    assert yc.rate_parameters(sub, Delta, normalize=False)["We"] == pytest.approx(rx["We"])
    with pytest.raises(ValueError):
        yc.rate_parameters(dict(omega=[0.5], weight=[1.0]), Delta)


def test_detailed_balance_degeneracy_factor():
    assert yc.detailed_balance(1e-3, 0.3, 0.0, g_ratio=2.0) == 0.0
    assert yc.detailed_balance(1e-3, 0.3, 0.1, g_ratio=2.0) == pytest.approx(2e-3*np.exp(-3.0))
