"""Regression tests for the `thermal` cluster of the 2026-09-25b hole hunt
(docs/audit_2026_09_25b_hole_hunt.md, findings 17 to 19), all in
thermal.anneal(), the imaginary-time stepper behind
Thermal_Spin_Chain.get_gs().

  17. anneal() returned early once one step changed the state by less than
      1e-7, a test of dtau^2*Var(H)/2 in absolute units: a Hamiltonian
      written below about 7e-3 came back as the T=infinity purification
      (E_th/s 0.000000 against -0.769460), and at unit scale a manifold
      split by 1e-2 was left unresolved (<Sz_tot> -0.089480 against
      -0.231059 at B=T=1e-2). The T>1e-5 switch to the plain ground state
      was absolute too.
  18. The steps were first-order, (1 - 0.1*H), on the unshifted H, so the
      error was set by 0.1 times the absolute, extensive energy: +12.3 per
      cent in E_th at n=10 and T=0.5, and at a +20 offset every factor was
      negative and the result sat above the T=infinity mean.
  19. nst = int(beta_half/0.1) is 0 above T=5: a float T raised
      ZeroDivisionError and a numpy scalar T returned the T=infinity state
      with a RuntimeWarning.

The stepper is now an order-k Taylor polynomial of exp(-dtau*(H-E_ref)),
E_ref the middle of the spectrum and dtau*W at most a fixed dimensionless
step, with no early return. So what is pinned is what a Gibbs state owes:
scale covariance E(sH, sT) = s E(H, T), shift covariance
E(H+c) - c = E(H), agreement with the exact Boltzmann average at a size-
independent tolerance, a convergence order, and n >= 1 steps at every T.
Anchors are exact Boltzmann averages over ED spectra of the physical chain.
"""

import io
import contextlib
import warnings
from math import factorial

import numpy as np
import pytest

from dmrgpy import cppext, spinchain, thermal


def _backend(version):
    marks = ()
    if version in (2, 3):
        marks = pytest.mark.skipif(not cppext.available(version),
            reason="needs the compiled itensor_version=%d extension" % version)
    ident = "v%d" % version if version in (2, 3) else str(version)
    return pytest.param(version, id=ident, marks=marks)


def _heis(ch, n, s=1.0, B=0.0, c=0.0):
    h = 0
    for i in range(n-1):
        h = h + s*(ch.Sx[i]*ch.Sx[i+1] + ch.Sy[i]*ch.Sy[i+1] + ch.Sz[i]*ch.Sz[i+1])
    if B != 0.0:
        for i in range(n): h = h + B*ch.Sz[i]
    if c != 0.0: h = h + c
    return h


def _spectrum(n, B=0.0):
    """ED spectrum of the physical chain, and <Sz_tot> in its eigenbasis"""
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version="python")
    sc.set_hamiltonian(_heis(sc, n, B=B))
    ed = sc.get_ED_obj()
    H = np.asarray(ed.get_hamiltonian().todense())
    sz = 0
    for i in range(n): sz = sz + sc.Sz[i]
    w, U = np.linalg.eigh(H)
    SZ = np.asarray(ed.MO2matrix(sz).todense())
    return w, np.real(np.diag(U.conj().T @ SZ @ U))


def _boltzmann(ev, T, obs=None):
    p = np.exp(-(ev - ev.min())/T)
    return np.sum((ev if obs is None else obs)*p)/np.sum(p)


def _old_euler(ev, T):
    """The pre-fix stepper in closed form: int(beta_half/0.1) steps of
    (1 - dtau*H) on the unshifted H (the early return never fires on the
    unit-scale chains this is used on)"""
    bh = 1./(2.*T); nst = int(bh/0.1); dt = bh/nst
    lw = 2*nst*np.log(np.abs(1. - dt*ev)); w = np.exp(lw - lw.max())
    return np.sum(ev*w)/np.sum(w)


def _new_stepper(ev, T, order, step):
    """What anneal() computes, in closed form over the exact spectrum"""
    bh = 1./(2.*T); W = ev.max() - ev.min()
    nst = max(1, int(np.ceil(bh*W/step - 1e-9))); dt = bh/nst
    x = dt*(ev - 0.5*(ev.max() + ev.min()))
    p = sum((-x)**j/factorial(j) for j in range(order+1))
    lw = 2*nst*np.log(np.abs(p)); w = np.exp(lw - lw.max())
    return np.sum(ev*w)/np.sum(w)


def _thermal(n, T, version="python", mode="ED", s=1.0, B=0.0, c=0.0,
             maxm=40, **attrs):
    """get_gs() of the purified chain; returns (energy - c, <Sz_tot>,
    number of annealing steps printed, warnings raised)"""
    np.random.seed(3)
    tc = thermal.Thermal_Spin_Chain(["S=1/2"]*n, itensor_version=version,
                                    mode=mode)
    tc.MBChain.maxm, tc.MBChain.nsweeps = maxm, 10
    for k, v in attrs.items(): setattr(tc, k, v)
    h = _heis(tc, n, s=s, B=B, c=c)
    tc.set_hamiltonian(h)
    tc.T = T
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), warnings.catch_warnings(record=True) as wl:
        warnings.simplefilter("always")
        wf = tc.get_gs()
    sz = 0
    for i in range(n): sz = sz + tc.Sz[i]
    e = float(np.real(wf.dot(h*wf)/wf.dot(wf))) - c
    m = float(np.real(wf.dot(sz*wf)/wf.dot(wf)))
    nst = buf.getvalue().count("Annealing, energy")
    return e, m, nst, [w for w in wl if issubclass(w.category, RuntimeWarning)]


EXACT3 = -0.769460  # 3-site Heisenberg chain, T = 0.5 J


@pytest.fixture(scope="module")
def ev3():
    return _spectrum(3)[0]


# ------------------------------------------------------------------ 17

@pytest.mark.parametrize("mode", ["ED", "DMRG"])
@pytest.mark.parametrize("s", [1e-6, 1e-3, 1.0, 1e2])
def test_thermal_energy_is_scale_covariant(ev3, mode, s):
    """E_th(sH, T=0.5 s)/s was -0.754223 at s=1, the T=infinity 0.000000 at
    s=1e-3 (the early return fired at the first step) and the ground
    state's -1.000000 at s=1e-6 (T=5e-7 fell under the absolute T>1e-5
    switch); the same state at every s, and the Boltzmann value"""
    e, m, nst, wl = _thermal(3, 0.5*s, mode=mode, s=s)
    assert _boltzmann(ev3, 0.5) == pytest.approx(EXACT3, abs=1e-6)
    assert e/s == pytest.approx(_boltzmann(ev3, 0.5), abs=1e-6)
    assert nst >= 1


@pytest.mark.parametrize("mode", ["ED", "DMRG"])
def test_a_manifold_split_by_the_temperature_is_resolved(mode):
    """3-site chain in a field B at T=B=1e-2: the doublet split by B is the
    whole calculation, and one step there changes the state by less than
    1e-7, so the early return fired after 200 of 500 steps (<Sz_tot>
    -0.089480 against -0.231059)"""
    B = 1e-2
    ev, szd = _spectrum(3, B=B)
    e, m, nst, wl = _thermal(3, B, mode=mode, B=B)
    assert _boltzmann(ev, B, szd) == pytest.approx(-0.231059, abs=1e-6)
    assert m == pytest.approx(_boltzmann(ev, B, szd), abs=2e-5)
    assert e == pytest.approx(_boltzmann(ev, B), abs=1e-6)


# ------------------------------------------------------------------ 18

@pytest.mark.parametrize("c", [-10.0, 20.0])
def test_thermal_energy_is_shift_covariant(ev3, c):
    """<H+c>-c was -0.422051 at c=-10 and +0.393848 at c=20 (above the
    T=infinity mean 0) against -0.769460 at every c; the Gibbs state does
    not depend on a constant"""
    assert _old_euler(ev3, 0.5) == pytest.approx(-0.754223, abs=1e-6)  # pre-fix
    e0 = _thermal(3, 0.5)[0]
    e, m, nst, wl = _thermal(3, 0.5, c=c)
    assert e == pytest.approx(e0, abs=1e-9)
    assert e == pytest.approx(_boltzmann(ev3, 0.5), abs=1e-6)


def test_shift_covariance_on_a_dmrg_backend(ev3):
    e, m, nst, wl = _thermal(3, 0.5, mode="DMRG", c=20.0)
    assert e == pytest.approx(_boltzmann(ev3, 0.5), abs=1e-6)


@pytest.mark.parametrize("n", [3, 4, 5, 6])
def test_error_does_not_grow_with_the_chain(n):
    """Relative error of E_th at T=0.5: the old stepper's closed form grows
    with n (+2.0, +4.2, +5.5, +7.0 per cent at n=3..6, kept here as the
    pre-fix reference), the new one stays below 1e-6"""
    ev = _spectrum(n)[0]
    exact = _boltzmann(ev, 0.5)
    old = abs(_old_euler(ev, 0.5) - exact)/abs(exact)
    assert old > 0.015*(n - 2)  # the defect, measured on the exact spectrum
    e, m, nst, wl = _thermal(n, 0.5)
    assert abs(e - exact)/abs(exact) < 1e-6


def test_error_does_not_grow_with_the_chain_on_dmrg():
    """n=8 on "python" DMRG, where the old stepper was 9.7 per cent off"""
    ev = _spectrum(8)[0]
    exact = _boltzmann(ev, 0.5)
    e, m, nst, wl = _thermal(8, 0.5, mode="DMRG", maxm=64)
    assert abs(e - exact)/abs(exact) < 1e-5


def test_get_gs_is_the_stepper_it_documents(ev3):
    """The implementation, independently of physics: on mode="ED" get_gs()
    is the closed form of an order-k Taylor step over the exact spectrum,
    with E_ref the middle of the band and nst = ceil(beta_half*W/step)"""
    ev = _spectrum(4)[0]
    for order, step in [(8, 2.0), (2, 0.5), (1, 0.25)]:
        e = _thermal(4, 0.5, anneal_order=order, anneal_step=step)[0]
        assert e == pytest.approx(_new_stepper(ev, 0.5, order, step), abs=1e-9)


@pytest.mark.parametrize("order", [1, 2, 4])
def test_the_error_converges_at_the_order_of_the_step(order):
    """Halving the dimensionless step divides the error by about 2^k (by
    (nst2/nst1)^k exactly, since the count is rounded up to an integer)"""
    ev = _spectrum(4)[0]
    exact = _boltzmann(ev, 0.5)
    step = {1: 0.2, 2: 0.4, 4: 0.5}[order]
    e1, _, n1, _ = _thermal(4, 0.5, anneal_order=order, anneal_step=step)
    e2, _, n2, _ = _thermal(4, 0.5, anneal_order=order, anneal_step=step/2)
    assert n2 == pytest.approx(2*n1, abs=1)
    assert abs((e1 - exact)/(e2 - exact)) == pytest.approx((n2/n1)**order, rel=0.25)


# ------------------------------------------------------------------ 19

@pytest.mark.parametrize("mode", ["ED", "DMRG"])
@pytest.mark.parametrize("T", [5.5, np.float64(7.0), np.linspace(1., 10., 4)[3],
                               np.float32(6.0), 50.0])
def test_every_finite_temperature_takes_a_step(ev3, mode, T):
    """T>5 took int(beta_half/0.1) = 0 steps: a float T raised
    ZeroDivisionError, a numpy scalar returned the T=infinity state
    (E_th 0.000000 against -0.055408 at T=7) with only a RuntimeWarning"""
    e, m, nst, wl = _thermal(3, T, mode=mode)
    assert nst >= 1
    assert not wl, [str(w.message) for w in wl]
    assert e == pytest.approx(_boltzmann(ev3, float(T)), abs=1e-7)


@pytest.mark.parametrize("J", [11.0, 100.0])
def test_large_couplings_at_half_their_scale(ev3, J):
    """T = J/2 raised ZeroDivisionError for every J > 10"""
    e = _thermal(3, 0.5*J, s=J)[0]
    assert e/J == pytest.approx(_boltzmann(ev3, 0.5), abs=1e-6)


def test_infinite_temperature_is_the_singlet_purification():
    e, m, nst, wl = _thermal(3, float("inf"))
    assert e == pytest.approx(0.0, abs=1e-12)  # Tr H / dim
    assert not wl


def test_a_temperature_set_after_construction_is_still_checked():
    tc = thermal.Thermal_Spin_Chain(["S=1/2"]*2, itensor_version="python",
                                    mode="ED")
    tc.set_hamiltonian(tc.Sz[0]*tc.Sz[1])
    for T in (-1.0, float("nan")):
        tc.T = T
        with pytest.raises(ValueError):
            tc.get_gs()


# ------------------------------------------------------------- backends

@pytest.mark.parametrize("version", [_backend(2), _backend(3), _backend("python")])
def test_every_backend_agrees_with_boltzmann(version):
    """n=4 is represented exactly at maxm=40, so what is left is the
    stepper, and that is size- and backend-independent"""
    ev = _spectrum(4)[0]
    for T in (0.5, 2.0):
        e = _thermal(4, T, version=version, mode="DMRG")[0]
        assert e == pytest.approx(_boltzmann(ev, T), abs=2e-6)
