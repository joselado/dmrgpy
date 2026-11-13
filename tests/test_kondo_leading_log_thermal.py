"""The leading-log Kondo spectrum at T>0 (leadinglog.leading_log_dIdV with
a thermalkernel.ThermalKernel). Energies in units of an exchange J, kB=1,
omega0=20 and Gamma0=2.5e-4, the values the construction was developed
with. Every reference is independent of the kernel tables: the thermal log
F(x,T) is a direct quadrature of w(u) F0(x+u), the third order is
conductance.py with its own FBuilder, and the T=0 route is the existing
flow. See docs/new_methods/kondo_leading_log/README.md, "Finite
temperature"."""
import warnings

import numpy as np
import pytest
from scipy.integrate import quad

from dmrgpy.kondospectrumtk import conductance, leadinglog
from dmrgpy.kondospectrumtk.stepfunctions import F0, FBuilder, Theta_prime
from dmrgpy.kondospectrumtk.thermalkernel import ThermalKernel

OM, G0 = 20., 2.5e-4


class KS():
    """A spectrum object from a dense Hamiltonian and the spin operators of
    the probed spin, with Boltzmann populations at T (kB=1)."""
    def __init__(self, H, S, T):
        e, v = np.linalg.eigh(H)
        self.e = e - e[0]
        self.dim = len(e)
        self.kB = 1.
        self.Sx, self.Sy, self.Sz = [(v.conj().T @ s @ v).astype(complex) for s in S]
        self.T = T
        p = (self.e <= 1e-9).astype(float) if T == 0. else np.exp(-self.e/T)
        self.p = p/p.sum()


def spin_ops(s):
    m = np.arange(s, -s - 1, -1)
    d = len(m)
    sp = np.zeros((d, d), dtype=complex)
    for i in range(1, d):
        sp[i-1, i] = np.sqrt(s*(s+1) - m[i]*(m[i]+1))
    return [(sp + sp.conj().T)/2, (sp - sp.conj().T)/2j, np.diag(m).astype(complex)]


def half(B, T):
    S = spin_ops(0.5)
    return KS(B*S[2], S, T)


def spin1(T):
    """S=1 with easy-axis and transverse anisotropy and a tilted field."""
    S = spin_ops(1.0)
    H = -1.0*S[2] @ S[2] + 0.2*(S[0] @ S[0] - S[1] @ S[1]) + 0.3*S[2] + 0.25*S[0]
    return KS(H, S, T)


def F_T(x, kT):
    """The thermal log F(x,T) = int du w(u) F0(x+u), by adaptive
    quadrature, independent of every table."""
    f = lambda u: Theta_prime(np.array([u/kT]))[0]/kT*F0(np.array([x + u]),
                                                        omega0=OM, Gamma0=G0)[0]
    pts = sorted({-x, -x - G0, -x + G0, -x - kT, -x + kT})
    pts = [p for p in pts if -60*kT < p < 60*kT]
    return quad(f, -60*kT, 60*kT, points=pts, limit=800, epsabs=1e-13,
                epsrel=1e-12)[0]


def ll(ks, eV, g, **kw):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return leadinglog.leading_log_dIdV(ks, np.asarray(eV, dtype=float), g,
                                           omega0=OM, Gamma0=G0, **kw)


def third(ks, eV, g, U, Fb):
    eV = np.asarray(eV, dtype=float)
    d = conductance.second_order_dIdV(ks, eV, U=U)
    d = d + conductance.third_order_kondo_dIdV(ks, eV, g, omega0=OM, Gamma0=G0, Fb=Fb)
    return d + conductance.third_order_potential_dIdV(ks, eV, g, U, omega0=OM,
                                                      Gamma0=G0, Fb=Fb)


@pytest.mark.parametrize("kT", [0.025, 0.3])
def test_the_kernel_integrates_to_the_thermal_log(kT):
    """int_0^inf dD K_T(x,D) = F(x,T): through the closed-form integral of
    the tables, and through K_T itself integrated over ln D."""
    tk = ThermalKernel(kT, OM, G0, 10*OM)
    for x in (0., G0, kT, 10*kT):
        ref = F_T(x, kT)
        assert abs(tk(np.array([x]), 1e10*OM, 0.)[0] - ref) < 1e-8*abs(ref)
    for x in (0., kT):
        I = quad(lambda s: np.exp(s)*tk(np.array([x]), np.exp(s))[0],
                 np.log(1e-10*G0), np.log(1e10*OM), limit=400, epsabs=1e-12,
                 epsrel=1e-11)[0]
        assert abs(I - F_T(x, kT)) < 1e-8*abs(I)


@pytest.mark.parametrize("kT", [0.025, 0.3])
def test_the_far_table_holds_where_the_thermal_part_is_not_negligible(kT):
    """Beyond xnear = 48 kT the kernel is K (1 + a spline in ln|x|); at these
    temperatures the thermal part there is 2e-5 to 1e-3 of F, so an error
    of the far table shows (measured below 8e-10)."""
    tk = ThermalKernel(kT, OM, G0, 10*OM)
    for c in (47.9, 48.1, 60., 150., 500.):
        x = c*kT
        ref = F_T(x, kT)
        assert abs(ref - F0(np.array([x]), OM, G0)[0]) > 1e-5*abs(ref)
        assert abs(tk(np.array([x]), 1e10*OM, 0.)[0] - ref) < 1e-8*abs(ref)


@pytest.mark.parametrize("model", [lambda T: half(1.0, T), spin1])
def test_the_zero_temperature_limit(model):
    """As kT/Gamma0 -> 0 the thermal route goes to the T=0 flow: as
    (kT/Gamma0)^2 over the first decade, more slowly below it."""
    eV = np.array([-2.3, -1.4, -0.6, -0.11, -0.012, 0.004, 0.03, 0.2, 0.55,
                   0.8, 1.7, 2.9])
    ref = ll(model(0.), eV, -0.05, U=0.25)
    got = ll(model(1e-4*G0), eV, -0.05, U=0.25)
    assert np.max(np.abs(got/ref - 1)) < 1e-9


@pytest.mark.parametrize("model,onshell", [(lambda T: half(1.0, T), "thermal"),
                                           (spin1, "thermal"),
                                           (lambda T: half(1.0, T), "T0")])
def test_the_first_order_is_the_thermal_third_order(model, onshell):
    """D(g) = LL - (second + third order at the same T). Halving g divides
    its even part by 4 and its odd part by 8; a mismatch at first order
    would give 2 in the odd part. The on-shell kernel enters at second
    order only, so onshell="T0" passes too."""
    kT, U = 0.3, 0.25
    ks = model(kT)
    eV = np.array([-2.5, -1.0, -0.4, 0.05, 0.7, 1.0, 1.3, 3.0])
    Fb = FBuilder(kT, omega0=OM, Gamma0=G0, kB=1.0, npts=32001)
    D = {g: ll(ks, eV, g, U=U, onshell=onshell) - third(ks, eV, g, U, Fb)
         for g in (0.01, -0.01, 0.005, -0.005)}
    even = lambda g: np.linalg.norm(D[g] + D[-g])/2
    odd = lambda g: np.linalg.norm(D[g] - D[-g])/2
    assert abs(even(0.01)/even(0.005) - 4) < 0.03
    assert abs(odd(0.01)/odd(0.005) - 8) < 0.08


@pytest.mark.parametrize("kT", [0.025, 0.3])
def test_a_free_spin_is_the_closed_form(kT):
    """A spin 1/2 at B=0: every process is elastic and one coupling flows,
    and with the thermal kernel everywhere the zero-bias peak is
    (3/4)/(1+g F(0,T))^2 in units of 2 pi."""
    g = -0.05
    d = ll(half(0., kT), [0.], g)[0]/(2*np.pi)
    closed = 0.75/(1 + g*F_T(0., kT))**2
    assert abs(d/closed - 1) < 1e-8


def test_the_stop_at_the_kondo_scale_is_the_zero_temperature_one():
    """A free spin at Jrho_s=-0.12 stops at the Kondo scale; at kT far below
    Gamma0 the stop, the NaN biases and the finite values are those of
    T=0."""
    eV = np.array([-3., -1., -0.3, -0.05, 0.001, 0.02, 0.2, 0.7, 2.])
    ref, iref = ll(half(0., 0.), eV, -0.12, U=0.25, return_info=True)
    got, info = ll(half(0., 1e-4*G0), eV, -0.12, U=0.25, return_info=True)
    assert iref["D_stop"] is not None
    assert abs(info["D_stop"]/iref["D_stop"] - 1) < 1e-6
    ok = ~np.isnan(ref)
    assert np.array_equal(ok, ~np.isnan(got))
    assert ok.any() and (~ok).any()
    assert np.max(np.abs(got[ok]/ref[ok] - 1)) < 1e-7


def test_bias_chunks_give_the_one_flow_result():
    """mem_gb cuts the biases in chunks, each with its own flow; the result
    is that of one flow, and a budget too small for 8 biases raises."""
    ks = spin1(0.3)
    eV = np.linspace(-2.4, 2.4, 24) + 1e-3
    one, i1 = ll(ks, eV, -0.05, U=0.25, return_info=True)
    # 8 biases per flow by the count of the route the flow takes, so 3 flows
    real = leadinglog._block_is_real(leadinglog._spin_ops(ks))
    fixed, per_w = leadinglog.memory_model(len(ks.e), len(conductance._occupied_states(ks)), real)
    cut, i3 = ll(ks, eV, -0.05, U=0.25, return_info=True,
                 mem_gb=(fixed + 16*per_w)*1.01/1e9)
    assert i1["n_chunks"] == 1 and i3["n_chunks"] == 3
    assert np.max(np.abs(cut/one - 1)) < 1e-10
    with pytest.raises(ValueError):
        ll(ks, eV, -0.05, mem_gb=1e-5)


def test_a_field_along_y_takes_the_complex_route_and_gives_the_same_curve():
    """A spin 1/2 in a field along y has a complex eigenbasis (the complex
    route); the tip is isotropic, so the curve is the one of the same field
    along z (the real route), up to the integrator's step sequence."""
    S = spin_ops(0.5)
    eV = np.array([-2.3, -0.6, -0.012, 0.03, 0.55, 1.7])
    kz, ky = KS(1.0*S[2], S, 0.3), KS(1.0*S[1], S, 0.3)
    assert not leadinglog._block_is_real(leadinglog._spin_ops(ky))
    a, b = ll(kz, eV, -0.05, U=0.25), ll(ky, eV, -0.05, U=0.25)
    assert np.max(np.abs(b/a - 1)) < 1e-9
