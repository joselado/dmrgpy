"""The leading-log resummation of the STM Kondo spectrum
(kondospectrumtk/leadinglog.py, docs/new_methods/kondo_leading_log/).

Three checks, in the order of the method's validation ladder: the closed
form at S=1/2, B=0; the S=1/2 equations in a field of Rosch, Costi,
Paaske, Woelfle, cond-mat/0301106 eq. 6, integrated independently; and
the first order in Jrho_s against conductance.py's third order on the
models of test_kondo_spectrum_tmatrix.py. Plus the stop at the Kondo
scale.
"""
import warnings
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.integrate import quad, solve_ivp

from dmrgpy.kondospectrumtk import conductance, leadinglog
from test_kondo_spectrum_tmatrix import MODELS, EVS, _ks, _spin

W0, G0 = 20e-3, 5e-6


def _K(x, D):
    return leadinglog._kernel(x, D, W0, G0)


def test_spin_half_at_zero_field_has_the_closed_form():
    """Every process is elastic, so the sample coupling runs as
    g(D) = g/(1 + g L(D)), L(D) = ln((D+omega0)/(D+Gamma0)), the exchange
    part of the tip vertex as t = exp(-int g(D) K(eV,D) dD), U does not
    run, and dI/dV = 2 pi (3/4 t^2 + 4 U^2)."""
    S = _spin(0.5)
    ks = SimpleNamespace(e=np.zeros(2), p=np.array([.5, .5]), dim=2, T=0.,
                         kB=8.617333262e-5, Sx=S[0], Sy=S[1], Sz=S[2])
    g, U = -0.05, 0.25
    eVs = np.array([-3e-3, -2e-4, 1e-5, 4e-4, 2e-3, 9e-3])
    got = leadinglog.leading_log_dIdV(ks, eVs, g, U=U, omega0=W0, Gamma0=G0)

    def gD(D): return g/(1. + g*np.log((D + W0)/(D + G0)))

    def t(v):
        f = lambda s: gD(np.exp(s))*_K(v, np.exp(s))*np.exp(s)
        val = quad(f, np.log(1e-12), np.log(1e6), limit=400, epsabs=1e-13,
                   points=[np.log(abs(v) + 1e-12), np.log(W0)])[0]
        return np.exp(-val)
    ref = np.array([2*np.pi*(0.75*t(v)**2 + 4*U**2) for v in eVs])
    got2, info = leadinglog.leading_log_dIdV(ks, eVs, g, U=U, omega0=W0,
                                             Gamma0=G0, return_info=True)
    assert info["D_stop"] is None
    # the coupling ran from |g| to |g/(1+g ln(omega0/Gamma0))|
    assert info["max_coupling"] == pytest.approx(abs(g/(1 + g*np.log(W0/G0))), rel=1e-6)
    assert np.max(np.abs(got - ref)) < 1e-7*np.max(ref)
    # the Kondo sign: Jrho_s < 0 makes the tip vertex grow
    assert np.all(got > 2*np.pi*(0.75 + 4*U**2))


def test_spin_half_in_a_field_follows_the_published_equations():
    """cond-mat/0301106 eq. 6 at omega = 0 and +-B/2 is closed in g_z(0)
    and g_perp(B/2): dg_z/dlnD = -2 g_perp^2 Theta_B and
    dg_perp/dlnD = -g_perp g_z (Theta_B + Theta_0), with
    g_theirs = -g_ours/2 (their tau is the Pauli matrix and their AF sign
    is positive). With the sharp Theta_x replaced by the same smooth
    kernel, D*K(x,D), the flowed sample vertex must follow it."""
    B = 1e-3
    S = _spin(0.5)
    ks = SimpleNamespace(e=np.array([0., B]), p=np.array([1., 0.]), dim=2, T=0.,
                         kB=8.617333262e-5, Sx=S[0], Sy=S[1], Sz=S[2])
    # eigenbasis of H = B Sz: state 0 is spin down, state 1 spin up
    P = np.array([[0., 1.], [1., 0.]])
    ks.Sx, ks.Sy, ks.Sz = [P @ x @ P for x in S]
    ks.e = np.array([0., B])
    g = -0.08
    Dmax, Dmin = 1e10*(W0 + B), 1e-10*G0
    sol, D_stop = leadinglog.sample_flow(ks, g, omega0=W0, Gamma0=G0,
                                         Dmax=Dmax, Dmin=Dmin)
    assert D_stop is None

    def rhs(s, y):
        D = np.exp(s); gz, gp = y
        return [-2*gp**2*D*_K(B, D), -gp*gz*D*(_K(B, D) + _K(0., D))]
    ref = solve_ivp(rhs, (np.log(Dmax), np.log(Dmin)), [-g/2, -g/2],
                    method="DOP853", rtol=1e-11, atol=1e-14, dense_output=True)
    for s in np.linspace(np.log(Dmax), np.log(Dmin), 40):
        G = sol.sol(s).view(complex).reshape(4, 2, 2)
        gz0 = G[3, 0, 0]/(ks.Sz[0, 0]/2); gz1 = G[3, 1, 1]/(ks.Sz[1, 1]/2)
        gpx = G[1, 0, 1]/(ks.Sx[0, 1]/2); gpy = G[2, 0, 1]/(ks.Sy[0, 1]/2)
        gz_t, gp_t = ref.sol(s)
        for x in (gz0, gz1):
            assert abs(-x/2 - gz_t) < 1e-8
        for x in (gpx, gpy):
            assert abs(-x/2 - gp_t) < 1e-8
        assert np.max(np.abs(G[0])) < 1e-10 # no potential scattering at omega=0


@pytest.mark.parametrize("name", list(MODELS))
def test_first_order_is_the_third_order(name):
    """Against second + third order of conductance.py, with U on: the
    difference is O(g^2), and its odd part in g is O(g^3), so the
    first order is exactly the third order (sign of g included)."""
    H, S = MODELS[name]
    ks = _ks(H, S)
    U = 0.25

    def LL(g): return leadinglog.leading_log_dIdV(ks, EVS, g, U=U, omega0=W0, Gamma0=G0)

    def third(g):
        return (conductance.third_order_kondo_dIdV(ks, EVS, g, omega0=W0, Gamma0=G0)
                + conductance.third_order_potential_dIdV(ks, EVS, g, U, omega0=W0, Gamma0=G0))
    d2 = conductance.second_order_dIdV(ks, EVS, U=U)
    r = {g: LL(g) - d2 - third(g) for g in (-0.01, -0.005, 0.005, 0.01)}
    ratio = np.max(np.abs(r[-0.01]))/np.max(np.abs(r[-0.005]))
    assert 3.5 < ratio < 4.6 # O(g^2)
    odd1 = (r[0.01] - r[-0.01])/2; odd2 = (r[0.005] - r[-0.005])/2
    assert np.max(np.abs(odd1)) < 0.05*np.max(np.abs(third(0.01)))
    assert 6. < np.max(np.abs(odd1))/np.max(np.abs(odd2)) < 10. # O(g^3)


def test_the_flow_stops_at_the_kondo_scale():
    """S=1/2 at B=0 with T_K above Gamma0: the coupling reaches order one,
    the biases near zero are NaN with a warning, the far ones are finite."""
    S = _spin(0.5)
    ks = SimpleNamespace(e=np.zeros(2), p=np.array([.5, .5]), dim=2, T=0.,
                         kB=8.617333262e-5, Sx=S[0], Sy=S[1], Sz=S[2])
    eVs = np.array([-1e-2, -1e-4, 0., 1e-4, 1e-2])
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        got, info = leadinglog.leading_log_dIdV(ks, eVs, -0.2, omega0=W0,
                                                Gamma0=G0, return_info=True)
    assert info["D_stop"] is not None and info["D_stop"] > G0
    assert info["max_coupling"] >= 1.
    assert any("order one" in str(x.message) for x in w)
    assert np.all(np.isnan(got[1:4]))
    assert np.all(np.isfinite(got[[0, 4]]))


def test_finite_temperature_is_refused():
    S = _spin(0.5)
    ks = SimpleNamespace(e=np.zeros(2), p=np.array([.5, .5]), dim=2, T=1.,
                         kB=8.617333262e-5, Sx=S[0], Sy=S[1], Sz=S[2])
    with pytest.raises(NotImplementedError):
        leadinglog.leading_log_dIdV(ks, [1e-3], -0.05)


def test_get_kondo_spectrum_order_LL():
    """order="LL" on a 3-site S=1 chain with anisotropy and a tilted field
    is the module's result on the chain's KondoSpectrum, is close to
    order=3 at small Jrho_s, and refuses T>0 and mode="DMRG"."""
    from dmrgpy import spinchain
    from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
    mub = 5.7883818060e-5
    sc = spinchain.Spin_Chain(["1", "1/2", "1"])
    h = -1e-3*sc.Sz[0]*sc.Sz[0] + 2*mub*(sc.Sx[0] + 2*sc.Sz[0])
    h = h + 0.5e-3*(sc.Sx[0]*sc.Sx[1] + sc.Sy[0]*sc.Sy[1] + sc.Sz[0]*sc.Sz[1])
    h = h + 0.3e-3*(sc.Sx[1]*sc.Sx[2] + sc.Sy[1]*sc.Sy[2] + sc.Sz[1]*sc.Sz[2])
    sc.set_hamiltonian(h)
    eVs = np.linspace(-3e-3, 3e-3, 13) + 1.3e-5
    _, got = sc.get_kondo_spectrum(eVs, site=0, Jrho_s=-0.05, U=0.25, T=0.,
                                   order="LL")
    ks = KondoSpectrum(sc, 0, 0.)
    ref = leadinglog.leading_log_dIdV(ks, eVs, -0.05, U=0.25)
    assert np.max(np.abs(got - ref)) < 1e-12*np.max(np.abs(ref))
    _, ll = sc.get_kondo_spectrum(eVs, site=0, Jrho_s=-0.002, U=0.25, T=0., order="LL")
    _, o3 = sc.get_kondo_spectrum(eVs, site=0, Jrho_s=-0.002, U=0.25, T=0., order=3)
    _, o2 = sc.get_kondo_spectrum(eVs, site=0, Jrho_s=-0.002, U=0.25, T=0., order=2)
    assert np.max(np.abs(ll - o3)) < 0.05*np.max(np.abs(o3 - o2))
    with pytest.raises(NotImplementedError):
        sc.get_kondo_spectrum(eVs, site=0, Jrho_s=-0.05, T=1.0, order="LL")
    with pytest.raises(NotImplementedError):
        sc.get_kondo_spectrum(eVs, site=0, Jrho_s=-0.05, T=0., order="LL", mode="DMRG")
    with pytest.raises(ValueError):
        sc.get_kondo_spectrum(eVs, site=0, order=4)
