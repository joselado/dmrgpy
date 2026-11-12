"""Regressions of the 2026-11 Kondo hole hunt, finding task-route-2: at T>0
the TASK unconverged-pole warning weighs each loose pole by how far it can
still move, d = min(r, r^2/g), against the width of the thermal step. On a
gapped chain where nkry=60 and nkry=200 agree to about 2e-6 of the peak it
no longer warns of a 3e-3 change, at T=0 nothing changes, and the size stays
between the true change and the old size."""
import warnings

import numpy as np
import pytest

from dmrgpy import spinchain
from dmrgpy.kondospectrumtk import conductance
from dmrgpy.kondospectrumtk.taskkondo import (CONVERGED_RESIDUAL,
                                              THERMAL_WIDTH, UNCONVERGED_STEP,
                                              TaskKondoSpectrum)

EMAX = 6e-3
SITE = 3


@pytest.fixture(scope="module")
def gapped():
    """A disordered XXZ chain with a staggered field, in eV: a gapped
    spectrum whose poles are far apart compared with kT at 0.1 K."""
    np.random.seed(1)
    n = 12
    rng = np.random.default_rng(5)
    sc = spinchain.Spin_Chain(["S=1/2"] * n, itensor_version="python")
    J = 1e-3 * rng.uniform(0.5, 1.5, n - 1)
    h = 0
    for i in range(n - 1):
        h = h + J[i] * (sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1]
                        + 1.5 * sc.Sz[i] * sc.Sz[i + 1])
    for i in range(n):
        h = h + (6e-4 * (-1) ** i + 1e-4 * rng.uniform(-1, 1)) * sc.Sz[i]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps, sc.cutoff = 24, 10, 1e-10
    sc.gs_energy()
    return sc


def _spectrum(sc, nkry, T):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return TaskKondoSpectrum(sc, SITE, nkry=nkry, T=T)


def _old_step(ks, emax):
    """The "step" of convergence() before 2026-11-12, kept verbatim: the
    largest weight of a loose pole, with no estimate of how far it moves."""
    w = sum(np.abs(M[:, 0]) ** 2 for M in (ks.Sx, ks.Sy, ks.Sz))
    total = max(float(np.sum(w)), 1e-300)
    width = max(float(np.max(np.abs(ks.e))), 1e-300)
    window = (ks.e > 0.) & (ks.e <= emax) & (w > 1e-14 * total)
    loose = window & (ks.residual > CONVERGED_RESIDUAL * width)
    if ks.T > 0. and loose.any():
        ew = np.sort(ks.e[w > 1e-14 * total])
        k = min(3, len(ew) - 1)
        apart = np.array([np.sort(np.abs(ew - x))[k] if k > 0 else np.inf
                          for x in ks.e])
        loose &= apart > THERMAL_WIDTH * ks.kB * ks.T
    return float(np.max(w[loose])) / total if loose.any() else 0.


def _true_change(ks, ref):
    """The largest change of the second-order dI/dV against a longer Krylov
    space on the same state, at every pole of both and on a fine grid."""
    pts = np.concatenate([ks.e[(ks.e > 0) & (ks.e < EMAX)],
                          ref.e[(ref.e > 0) & (ref.e < EMAX)],
                          np.linspace(0., EMAX, 3001)])
    pts = np.unique(pts)
    a = conductance.second_order_dIdV(ks, pts)
    b = conductance.second_order_dIdV(ref, pts)
    return np.max(np.abs(a - b)) / np.max(np.abs(b))


@pytest.mark.parametrize("nkry", [20, 60])
def test_the_zero_temperature_step_is_unchanged(gapped, nkry):
    ks = _spectrum(gapped, nkry, 0.)
    assert ks.convergence(EMAX)["step"] == _old_step(ks, EMAX)


def test_a_gapped_chain_at_finite_temperature_no_longer_warns(gapped):
    ks = _spectrum(gapped, 60, 0.1)
    ref = _spectrum(gapped, 200, 0.1)
    true = _true_change(ks, ref)
    assert true < 1e-4  # nkry=60 is converged at the scale of the warning
    assert _old_step(ks, EMAX) > UNCONVERGED_STEP  # the old false alarm
    assert ks.convergence(EMAX)["step"] < UNCONVERGED_STEP
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        ks.warn_unconverged(EMAX)


@pytest.mark.parametrize("nkry", [20, 30, 40, 60])
def test_the_size_stays_between_the_true_change_and_the_old_size(gapped,
                                                                   nkry):
    ks = _spectrum(gapped, nkry, 0.1)
    ref = _spectrum(gapped, 200, 0.1)
    step = ks.convergence(EMAX)["step"]
    assert step <= _old_step(ks, EMAX)
    assert step >= _true_change(ks, ref)


def test_the_shift_is_quadratic_in_the_residual_and_linear_when_degenerate():
    ks = object.__new__(TaskKondoSpectrum)
    ks.e = np.array([0., 1e-3, 1e-3, 2e-3])
    ks.residual = np.array([0., 1e-5, 1e-5, 1e-5])
    ks.kB, ks.T = 1., 1e-5
    f = ks._thermal_shift_factor(np.array([1, 2, 3]))
    width = THERMAL_WIDTH * 1e-5
    # rows 1 and 2 are degenerate, so only d <= r is known
    assert f[0] == pytest.approx(1e-5 / width)
    assert f[1] == pytest.approx(1e-5 / width)
    # row 3 is 1e-3 from its neighbour: d = r^2/g
    assert f[2] == pytest.approx((1e-5) ** 2 / 1e-3 / width)
    ks.T = 1e-9  # a step far narrower than any shift: the full weight
    assert np.all(ks._thermal_shift_factor(np.array([1, 2, 3])) == 1.)
