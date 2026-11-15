"""Purification by implicit steps on a logarithmic grid
(Thermal_Spin_Chain.anneal_method = "implicit", thermal.anneal_implicit).

The reference is the exact Boltzmann average over the ED spectrum of the
physical chain, as in test_audit_2026_09_25b_thermal.py, whose helpers are
reused. The closed-form test pins the reason for the default (order=8,
nsteps=2): the purification starts from the T=infinity state, and the
trapezoid rule at the METTS default is too coarse there.
"""

import numpy as np
import pytest

from dmrgpy import cppext, imaginarytime, thermal

from test_audit_2026_09_25b_thermal import _spectrum, _boltzmann, _thermal
from _helpers import julia_live_param

# every julia_live test carries the julia_live id, so -k "not julia_live"
# skips them
JL = pytest.mark.parametrize("version", [julia_live_param()])


def _closed_form(ev, T, order, nsteps, tau0=0.01):
    """The thermal energy of the purified state after the rational
    propagator of the grid, over the exact spectrum"""
    x = ev - ev.min()
    lw = np.zeros_like(ev)
    for dt in imaginarytime.log_grid(1./(2.*T), tau0, nsteps):
        for b in imaginarytime.step_factors(dt, order):
            lw += 2*np.log(np.abs((1 - b*x)/(1 + b*x)))
    w = np.exp(lw - lw.max())
    return np.sum(ev*w)/np.sum(w)


def test_the_default_rule_is_accurate_from_infinite_temperature():
    """6-site chain at T=1, the worst temperature for the time step: the
    trapezoid rule with 8 steps per panel is 4.6e-4 off, the default
    order=8 with 2 steps per panel below 1e-8"""
    ev = _spectrum(6)[0]
    exact = _boltzmann(ev, 1.0)
    trap = abs(_closed_form(ev, 1.0, 2, 8) - exact)/abs(exact)
    gauss = abs(_closed_form(ev, 1.0, 8, 2) - exact)/abs(exact)
    assert trap > 1e-4
    assert gauss < 1e-8


@pytest.mark.parametrize("T", [5.0, 0.5, 0.05])
def test_ed_route_is_the_boltzmann_average(T):
    """mode="ED" evolves the exact e^{-tau H} on the grid"""
    ev = _spectrum(4)[0]
    e, m, nst, wl = _thermal(4, T, mode="ED", anneal_method="implicit")
    assert nst == 0  # no Taylor step was taken
    assert e == pytest.approx(_boltzmann(ev, T), abs=1e-9)


@pytest.mark.parametrize("version", [
    "python",
    pytest.param(3, marks=pytest.mark.skipif(not cppext.available(3),
                 reason="needs the compiled itensor_version=3 extension"))])
@pytest.mark.parametrize("T", [1.0, 0.1])
def test_dmrg_route_is_the_boltzmann_average(version, T):
    """4.1e-6 (relative) at T=1 and 6e-8 at T=0.1 on both backends: the
    floor is the truncation at cutoff=1e-12 of every factor solve, four
    per step at order=8, not the time step (below 1e-8 here)"""
    ev = _spectrum(4)[0]
    exact = _boltzmann(ev, T)
    e, m, nst, wl = _thermal(4, T, version=version, mode="DMRG",
                             anneal_method="implicit")
    assert nst == 0
    assert abs(e - exact)/abs(exact) < 1e-5


@JL
def test_julia_live_route_is_the_boltzmann_average(version):
    """the DMRG route on julia_live (mpsjulialive/implicit.jl) at T=0.1,
    within the bound of the "python" and v3 routes"""
    ev = _spectrum(4)[0]
    exact = _boltzmann(ev, 0.1)
    e, m, nst, wl = _thermal(4, 0.1, version=version, mode="DMRG",
                             anneal_method="implicit")
    assert nst == 0
    assert abs(e - exact)/abs(exact) < 1e-5


def test_options_reach_the_steps():
    """the trapezoid rule through anneal_options gives its own, larger
    error, the one of its closed form"""
    ev = _spectrum(4)[0]
    e, m, nst, wl = _thermal(4, 1.0, mode="DMRG", anneal_method="implicit",
                             anneal_options=dict(order=2, nsteps=8))
    assert e == pytest.approx(_closed_form(ev, 1.0, 2, 8), abs=1e-5)
    assert abs(e - _boltzmann(ev, 1.0)) > 1e-4


def test_infinite_temperature_is_the_singlet_purification():
    e, m, nst, wl = _thermal(3, float("inf"), anneal_method="implicit")
    assert e == pytest.approx(0.0, abs=1e-12)


def test_unknown_method_and_option_raise():
    with pytest.raises(ValueError):
        _thermal(3, 0.5, anneal_method="tdvp")
    with pytest.raises(ValueError):  # options the Taylor stepper ignores
        _thermal(3, 0.5, anneal_options=dict(order=4))
    with pytest.raises(TypeError):
        _thermal(3, 0.5, mode="DMRG", anneal_method="implicit",
                 anneal_options=dict(shift=0.0))


@pytest.mark.skipif(not cppext.available(2),
                    reason="needs the compiled itensor_version=2 extension")
def test_v2_raises():
    with pytest.raises(NotImplementedError):
        _thermal(4, 0.5, version=2, mode="DMRG", anneal_method="implicit")


# ---------------------------------------------------- the default, None

_v3 = pytest.mark.skipif(not cppext.available(3),
                         reason="needs the compiled itensor_version=3 extension")


def test_counts_behind_the_choice():
    """4 sites: 480 Taylor applications against 112 solves at T=0.01, the
    grid; 16 against 64 at T=0.5, Taylor"""
    tc = thermal.Thermal_Spin_Chain(["S=1/2"]*4, itensor_version="python")
    h = sum(tc.Sx[i]*tc.Sx[i+1] + tc.Sy[i]*tc.Sy[i+1] + tc.Sz[i]*tc.Sz[i+1]
            for i in range(3))
    tc.set_hamiltonian(h)
    ev = _spectrum(4)[0]
    erange = (ev.min(), ev.max())
    ta, im = thermal.anneal_counts(tc, 0.01, erange)
    assert ta >= thermal._AUTO_RATIO*im
    ta, im = thermal.anneal_counts(tc, 0.5, erange)
    assert ta < im


@_v3
def test_default_takes_the_grid_on_v3_at_low_temperature():
    ev = _spectrum(4)[0]
    e, m, nst, wl = _thermal(4, 0.01, version=3, mode="DMRG")
    assert nst == 0  # no Taylor step was printed
    assert abs(e - _boltzmann(ev, 0.01))/abs(_boltzmann(ev, 0.01)) < 1e-6
    e, m, nst, wl = _thermal(4, 0.5, version=3, mode="DMRG")
    assert nst > 0


@JL
def test_default_keeps_taylor_on_julia_live(version):
    """the automatic choice does not take the grid on julia_live, at a T
    where v3 takes it (test_default_takes_the_grid_on_v3_at_low_temperature)"""
    tc = thermal.Thermal_Spin_Chain(["S=1/2"]*4, itensor_version=version)
    h = sum(tc.Sx[i]*tc.Sx[i+1] + tc.Sy[i]*tc.Sy[i+1] + tc.Sz[i]*tc.Sz[i+1]
            for i in range(3))
    tc.set_hamiltonian(h)
    ev = _spectrum(4)[0]
    wf = tc.MBChain.random_mps()
    assert thermal.auto_anneal_method(tc, wf, 0.01,
                                      (ev.min(), ev.max())) == "taylor"


@pytest.mark.parametrize("mode", ["ED", "DMRG"])
def test_default_keeps_taylor_on_python_and_ed(mode):
    e, m, nst, wl = _thermal(4, 0.01, version="python", mode=mode)
    assert nst > 0


@_v3
def test_default_keeps_a_configured_taylor_stepper():
    e, m, nst, wl = _thermal(4, 0.01, version=3, mode="DMRG", anneal_order=8)
    assert nst > 0


@_v3
def test_default_keeps_taylor_above_the_temperature_bound(monkeypatch):
    """with any count ratio accepted, T=0.5 > 0.05 J still takes Taylor:
    the count alone takes the grid at T=0.5 on 100 sites, where it is less
    accurate and was never timed"""
    monkeypatch.setattr(thermal, "_AUTO_RATIO", 0.0)
    e, m, nst, wl = _thermal(4, 0.5, version=3, mode="DMRG")
    assert nst > 0
    e, m, nst, wl = _thermal(4, 0.04, version=3, mode="DMRG")
    assert nst == 0
