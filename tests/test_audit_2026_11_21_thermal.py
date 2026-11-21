"""Regression tests for the `thermal` cluster of the 2026-11-21 hole hunt
(docs/audit_2026_11_21_hole_hunt.md, finding 10).

  10. Thermal_Spin_Chain.get_gs() returned its stored purified state
      whenever computed_gs was set, and only set_hamiltonian() cleared it.
      A temperature set after the first solve was ignored: on a 3-site
      Heisenberg chain the same object gave -0.206813 at T=0.2, the energy
      at T=2.0, against the exact -0.991668, on mode="ED" and on "python".
      anneal_step, anneal_order, anneal_method, anneal_options, mode and the
      solver settings of MBChain were ignored the same way, and an invalid
      setting after a solve was not refused.

The stored state now carries the key it was purified under
(thermal.purification_key) and comes back only while that key holds. What
is pinned: a temperature sweep on one object agrees with fresh objects and
with the exact Boltzmann average; an unchanged call returns the stored
state without a solve; an edit in place of anneal_options, a new stepper, a
new mode and a new maxm or nsweeps on MBChain purify again; an invalid
setting after a solve raises and leaves the stored state as it was. The
anchors are exact Boltzmann averages and closed forms of the two steppers
over the ED spectrum of the physical chain.
"""

import io
import contextlib

import numpy as np
import pytest

from dmrgpy import cppext, imaginarytime, thermal

from test_audit_2026_09_25b_thermal import _spectrum, _new_stepper
from test_thermal_anneal_implicit import _closed_form

N = 3  # physical sites; the purified chain has 2N


def _heis(ch, n):
    h = 0
    for i in range(n-1):
        h = h + ch.Sx[i]*ch.Sx[i+1] + ch.Sy[i]*ch.Sy[i+1] + ch.Sz[i]*ch.Sz[i+1]
    return h


def _exact(ev, T):
    """The Boltzmann average of H over the spectrum ev, T=0 included"""
    if T == 0.: return float(ev.min())
    w = np.exp(-(ev - ev.min())/T)
    return float(np.sum(ev*w)/np.sum(w))


def _make(T, mode, n=N, version="python", maxm=16, **attrs):
    """A purified chain at T with the Heisenberg Hamiltonian, seeded, with
    its sweep schedule pinned and the attributes attrs set on it"""
    np.random.seed(3)
    tc = thermal.Thermal_Spin_Chain(["S=1/2"]*n, T=T, mode=mode,
                                    itensor_version=version)
    tc.MBChain.maxm, tc.MBChain.nsweeps = maxm, 10
    for k, v in attrs.items(): setattr(tc, k, v)
    h = _heis(tc, n)
    tc.set_hamiltonian(h)
    return tc, h


def _solve(tc, h):
    """get_gs(), the energy <H> of its state, and the number of Taylor
    steps it printed (0 when no Taylor step ran)"""
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        wf = tc.get_gs()
    e = float(np.real(wf.dot(h*wf)/wf.dot(wf)))
    return wf, e, buf.getvalue().count("Annealing, energy")


def _quiet(f):
    with contextlib.redirect_stdout(io.StringIO()):
        return f()


@pytest.mark.parametrize("mode", ["ED", "DMRG"])
def test_a_temperature_sweep_on_one_object_follows_the_temperature(mode):
    """T=2.0, then 0.2, back to 2.0, then the plain ground state at T=0, on
    one object: each the energy of a fresh object at that T and the exact
    Boltzmann value. Before, every T after the first gave the T=2.0 energy,
    0.785 off at T=0.2. A second call with nothing changed returns the
    stored state itself, without a solve"""
    ev = _spectrum(N)[0]
    assert abs(_exact(ev, 2.0) - _exact(ev, 0.2)) > 0.5  # what was missed
    tc, h = _make(2.0, mode)
    for T in (2.0, 0.2, 2.0, 0.0):
        tc.T = T
        wf, e, nst = _solve(tc, h)
        fresh, hf = _make(T, mode)
        assert e == pytest.approx(_solve(fresh, hf)[1], abs=1e-7)
        assert e == pytest.approx(_exact(ev, T), abs=1e-6)
        again, e2, nst2 = _solve(tc, h)
        assert again is wf and nst2 == 0


def test_a_new_stepper_after_a_solve_is_the_one_taken():
    """anneal_order and anneal_step set after a solve: the state is the
    closed form of the new stepper, not of the old one"""
    ev = _spectrum(N)[0]
    default = _new_stepper(ev, 0.5, 8, 2.0)
    coarse = _new_stepper(ev, 0.5, 2, 0.5)
    assert abs(default - coarse) > 1e-5
    tc, h = _make(0.5, "ED")
    assert _solve(tc, h)[1] == pytest.approx(default, abs=1e-9)
    tc.anneal_order, tc.anneal_step = 2, 0.5
    assert _solve(tc, h)[1] == pytest.approx(coarse, abs=1e-9)


def test_an_edit_in_place_of_anneal_options_is_read_on_dmrg():
    """the trapezoid rule with 2 steps per panel, then the same dict edited
    in place to 4: the second state is the closed form of 4 steps"""
    ev = _spectrum(N)[0]
    opts = dict(order=2, nsteps=2)
    tc, h = _make(1.0, "DMRG", anneal_method="implicit", anneal_options=opts)
    tau0 = 0.01/imaginarytime._unit(tc.hamiltonian)  # the default of the grid
    two = _closed_form(ev, 1.0, 2, 2, tau0=tau0)
    four = _closed_form(ev, 1.0, 2, 4, tau0=tau0)
    gap = abs(two - four)  # about 5e-4 by the error of the trapezoid rule
    assert gap > 1e-4
    wf, e, nst = _solve(tc, h)
    assert nst == 0
    assert e == pytest.approx(two, abs=0.1*gap)
    opts["nsteps"] = 4  # tc.anneal_options is this dict
    wf2, e2, nst2 = _solve(tc, h)
    assert wf2 is not wf
    assert e2 == pytest.approx(four, abs=0.1*gap)


def test_an_edit_in_place_of_anneal_options_is_read_on_ed():
    """an unknown keyword put into the dict after a solve raises, as on a
    fresh object, and once it is removed the state is purified again, with
    MBChain holding it rather than the singlet the failed solve left"""
    ev = _spectrum(N)[0]
    opts = {}
    tc, h = _make(1.0, "ED", anneal_method="implicit", anneal_options=opts)
    wf, e, nst = _solve(tc, h)
    assert e == pytest.approx(_exact(ev, 1.0), abs=1e-9)
    opts["shift"] = 0.0
    with pytest.raises(TypeError):
        _quiet(tc.get_gs)
    del opts["shift"]
    wf2, e2, nst2 = _solve(tc, h)
    assert e2 == pytest.approx(e, abs=1e-9)
    assert np.real(_quiet(lambda: tc.MBChain.vev(h))) == pytest.approx(e, abs=1e-9)
    opts["nsteps"] = 3  # a valid edit is a new key too
    assert _solve(tc, h)[0] is not wf2


@pytest.mark.parametrize("name,value", [("maxm", 12), ("nsweeps", 8)])
def test_a_new_solver_setting_of_the_chain_purifies_again(name, value):
    """a setting of MBChain changed after a solve: a new purification at
    that setting (Taylor steps printed again), with the Boltzmann energy"""
    ev = _spectrum(N)[0]
    tc, h = _make(0.5, "DMRG")
    wf, e, nst = _solve(tc, h)
    assert nst > 0
    setattr(tc.MBChain, name, value)
    wf2, e2, nst2 = _solve(tc, h)
    assert wf2 is not wf and nst2 > 0
    assert e2 == pytest.approx(_exact(ev, 0.5), abs=1e-6)


def test_a_new_mode_purifies_again():
    """tc.mode="DMRG" after an ED solve gives an MPS, not the ED state"""
    ev = _spectrum(N)[0]
    tc, h = _make(0.5, "ED")
    wf, e, nst = _solve(tc, h)
    assert wf.mode == "ED"
    tc.mode = "DMRG"
    wf2, e2, nst2 = _solve(tc, h)
    assert wf2.mode == "DMRG" and nst2 > 0
    assert e2 == pytest.approx(_exact(ev, 0.5), abs=1e-6)


def test_an_invalid_setting_after_a_solve_raises():
    """refused on a solved object as on a fresh one, and the stored state
    is kept: putting the valid setting back returns it without a solve"""
    tc, h = _make(0.5, "ED")
    wf = _solve(tc, h)[0]
    for name, bad in [("anneal_method", "not_a_method"), ("T", float("nan")),
                      ("T", -1.0), ("anneal_options", dict(order=4))]:
        good = getattr(tc, name)
        setattr(tc, name, bad)
        with pytest.raises(ValueError):
            _quiet(tc.get_gs)
        setattr(tc, name, good)
        again, e, nst = _solve(tc, h)
        assert again is wf and nst == 0


def test_the_key_sees_a_stepper_set_on_the_instance():
    """auto_anneal_method() keeps the Taylor stepper once anneal_order is
    set on the instance, at the default value too, so the key changes"""
    tc, h = _make(0.5, "ED")
    key = thermal.purification_key(tc, 0.5)
    assert thermal.purification_key(tc, 0.5) == key
    tc.anneal_order = thermal.Thermal_Spin_Chain.anneal_order
    assert thermal.purification_key(tc, 0.5) != key


@pytest.mark.skipif(not cppext.available(3),
                    reason="needs the compiled itensor_version=3 extension")
def test_anneal_order_set_after_a_solve_keeps_the_taylor_stepper_on_v3():
    """4 sites at T=0.01 on v3: the default takes the implicit grid (no
    Taylor step printed), and anneal_order=8, its default value, set after
    the solve takes the Taylor stepper, as on a fresh object"""
    tc, h = _make(0.01, "DMRG", n=4, version=3, maxm=40)
    assert _solve(tc, h)[2] == 0
    tc.anneal_order = 8
    assert _solve(tc, h)[2] > 0
