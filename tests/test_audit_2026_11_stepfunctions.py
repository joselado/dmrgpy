"""Kondo step functions and ED eigenbasis (2026-11 hole hunt, cluster
stepfunctions): Theta_prime keeps its digits for small |x|, FBuilder
evaluates F0 only outside its table, cached_fbuilder shares one object,
and the sparse eigenbasis spin operator equals the old dense formula."""
import numpy as np
import pytest
import mpmath as mp

from dmrgpy import spinchain, fermionchain
from dmrgpy.algebra import algebra
from dmrgpy.kondospectrumtk import stepfunctions as sf
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum


def _tp_exact(x):
    mp.mp.dps = 50
    x = mp.mpf(float(x))
    h = x/2
    return float((x*mp.cosh(h) - 2*mp.sinh(h))/(4*mp.sinh(h)**3))


@pytest.mark.parametrize("x", [1e-4, 3e-4, 1e-3, 1e-2, 0.049, 0.051, 0.3, 2.,
                               -1e-3, -0.04])
def test_theta_prime_keeps_digits(x):
    got = float(sf.Theta_prime(np.array([x]))[0])
    assert abs(got - _tp_exact(x)) < 1e-12


def test_theta_prime_unit_normalized():
    u = np.linspace(-60, 60, 120001)
    assert abs(np.trapezoid(sf.Theta_prime(u), u) - 1) < 1e-8


def _old_call(fb, x):
    x = np.atleast_1d(np.asarray(x, dtype=float))
    inside = (x >= fb._F_min) & (x <= fb._F_max)
    out = sf.F0(x, omega0=fb.omega0, Gamma0=fb.Gamma0)
    out[inside] = fb._F_spline(x[inside])
    return out


def test_fbuilder_call_equals_old_overwrite_form():
    fb = sf.FBuilder(1., omega0=2e-3, npts=401, grid_span=4.)
    lo, hi = fb._F_min, fb._F_max
    inner = np.linspace(lo, hi, 50)
    outer = np.array([lo - 1e-3, lo - 1e-2, hi + 1e-3, hi + 1e-1])
    for x in (inner, outer, np.concatenate([outer[:2], inner, outer[2:]]),
              0.003):
        assert np.array_equal(fb(x), _old_call(fb, x))


def test_cached_fbuilder_is_shared_and_keyed():
    a = sf.cached_fbuilder(1., omega0=2e-3, npts=401, grid_span=4.)
    b = sf.cached_fbuilder(1, omega0=2e-3, npts=401, grid_span=4.)
    assert a is b
    assert sf.cached_fbuilder(1., omega0=3e-3, npts=401, grid_span=4.) is not a
    assert sf.cached_fbuilder(1., omega0=2e-3, kB=1e-4, npts=401,
                              grid_span=4.) is not a


def _old_eigenbasis(ks, name):
    op = getattr(ks.chain, name)[ks.site]
    mat = np.array(algebra.todense(ks.chain.get_ED_obj().get_operator(op)),
                   dtype=complex)
    return np.conjugate(ks.vs.T)@mat@ks.vs


def _check(chain):
    ks = KondoSpectrum(chain, 0, 0.)
    for name in ("Sx", "Sy", "Sz"):
        d = getattr(ks, name) - _old_eigenbasis(ks, name)
        assert np.max(np.abs(d)) < 1e-12


def test_eigenbasis_spin_operator_spin_chain():
    sc = spinchain.Spin_Chain(["1/2"]*5)
    h = 0
    for i in range(4):
        h = h + sum(a[i]*a[i+1] for a in (sc.Sx, sc.Sy, sc.Sz))
    sc.set_hamiltonian(h + 0.1*sc.Sz[0] + 0.05*sc.Sx[2])
    _check(sc)


def test_eigenbasis_spin_operator_hubbard_sector():
    fc = fermionchain.Spinful_Fermionic_Chain_Native(
        3, itensor_version="python")
    h = 0
    for i in range(2):
        h = h - (fc.Cdagup[i]*fc.Cup[i+1] + fc.Cdagdn[i]*fc.Cdn[i+1])
    h = h + h.get_dagger()
    for i in range(3):
        h = h + 2.*(fc.Nup[i]-0.5)*(fc.Ndn[i]-0.5) + 0.1*fc.Sz[i]
    fc.set_hamiltonian(h)
    fc.set_conserved_sector(Nf=3)
    _check(fc)
