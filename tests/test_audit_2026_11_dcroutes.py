"""The es coverage check of second_order_dIdV_dc (2026-11 Kondo hunt,
correlator-twotime-2). The lower edge cuts the elastic line when es starts
at 0 and the ground state has a moment; the top is a necessary condition.
Values do not change: a covered grid gives the same numbers as before."""
import warnings

import numpy as np
import pytest

from dmrgpy import spinchain
from dmrgpy.kondospectrumtk.secondorder_dc import second_order_dIdV_dc

G, MUB = 2.0, 5.7883818066e-5


def _chain(B):
    sc = spinchain.Spin_Chain(["1/2"]*4)
    h = G*MUB*B*sc.Sz[0]
    for i in range(3):
        h = h + 1e-3*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
                      + sc.Sz[i]*sc.Sz[i+1])
    sc.set_hamiltonian(h)
    return sc


def _run(sc, lo, hi):
    es = np.linspace(lo, hi, int((hi-lo)/5e-6))
    eVs = np.linspace(-4e-3, 4e-3, 9)
    with warnings.catch_warnings(record=True) as W:
        warnings.simplefilter("always")
        r = second_order_dIdV_dc(sc, 0, eVs, mode="ED", submode="ED",
                                 delta=2e-5, es=es)
    return r, [w for w in W if issubclass(w.category, RuntimeWarning)]


def test_es_from_zero_warns_and_cuts_the_elastic_line():
    sc = _chain(10.)
    bad, w = _run(sc, 0.0, 1e-2)
    good, w2 = _run(sc, -1e-3, 1e-2)
    assert len(w) >= 1 and "lowest point" in str(w[0].message)
    assert w2 == []
    assert abs(bad[4]) < 0.5*abs(good[4])


def test_covered_grid_does_not_warn():
    _, w = _run(_chain(10.), -1e-3, 1e-2)
    assert w == []


def test_singlet_with_es_from_zero_does_not_warn():
    _, w = _run(_chain(0.), 0.0, 1e-2)
    assert w == []


def test_es_below_max_eV_warns():
    _, w = _run(_chain(10.), -1e-3, 2e-3)
    assert any("below max|eV|" in str(x.message) for x in w)
