"""Item 13 of docs/kondo_speedup_plan.md and the memory warning of the
leading-log flow on the CPU: memory_model counts real numbers on the real
route, and a flow with mem_gb=None that the model puts above what the
process can still take warns and names mem_gb."""

import warnings

import numpy as np
import pytest

from dmrgpy import spinchain
from dmrgpy.kondospectrumtk import conductance, edkondo, leadinglog


def _chain(n=3, field_y=False):
    sc = spinchain.Spin_Chain(["S=1/2"]*n)
    h = 0
    for i in range(n - 1):
        h = h + 1e-3*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
    for i in range(n):
        h = h + 0.3e-3*sc.Sz[i] + 0.1e-3*(sc.Sy[i] if field_y else sc.Sx[i])
    sc.set_hamiltonian(h)
    return sc


def test_the_real_route_counts_half_the_bytes():
    for d, m in ((16, 1), (512, 4), (256, 223)):
        c = leadinglog.memory_model(d, m)
        r = leadinglog.memory_model(d, m, real=True)
        assert r[0]*2 == c[0] and r[1]*2 == c[1]
    # the complex count is the one the docstring measured: 7.0 GB at 9 spins
    fixed, per_w = leadinglog.memory_model(512, 4)
    assert (fixed + 1642*per_w)/1e9 == pytest.approx(7.0, abs=0.05)


@pytest.mark.parametrize("field_y,real", [(False, True), (True, False)])
def test_mem_gb_cuts_by_the_count_of_the_route(field_y, real, monkeypatch):
    """A budget of 8 biases by the count of the route gives three flows for
    24 biases, on the real route (a field along x) and on the complex one
    (along y)."""
    sc = _chain(field_y=field_y)
    ks = edkondo.KondoSpectrum(sc, 0, 2.0)
    assert leadinglog._block_is_real(leadinglog._spin_ops(ks)) == real
    occ = conductance._occupied_states(ks)
    fixed, per_w = leadinglog.memory_model(len(ks.e), len(occ), real)
    eV = np.linspace(-3e-3, 3e-3, 24) + 1e-5
    _, info = leadinglog.leading_log_dIdV(ks, eV, -0.05, U=0.2, return_info=True,
                                          mem_gb=(fixed + 16*per_w)*1.01/1e9)
    assert info["n_chunks"] == 3


def test_a_flow_above_the_usable_memory_warns(monkeypatch):
    sc = _chain()
    ks = edkondo.KondoSpectrum(sc, 0, 0.)
    eV = np.linspace(-3e-3, 3e-3, 11) + 1e-5
    monkeypatch.setattr(leadinglog, "usable_memory", lambda: 1e3)
    with pytest.warns(RuntimeWarning, match="mem_gb"):
        ref = leadinglog.leading_log_dIdV(ks, eV, -0.05, U=0.2)
    for usable, kw in ((lambda: 1e3, dict(mem_gb=1.)), (lambda: None, {}),
                       (lambda: 1e15, {})):
        monkeypatch.setattr(leadinglog, "usable_memory", usable)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            got = leadinglog.leading_log_dIdV(ks, eV, -0.05, U=0.2, **kw)
        assert not any("can still take" in str(x.message) for x in w)
        assert np.max(np.abs(got - ref)) <= 1e-12*np.max(np.abs(ref))


def test_usable_memory_reads_this_process():
    free = leadinglog.usable_memory()
    assert free is None or free > 0
