"""Coverage for the T=0 two-time-correlator construction of the
third-order Kondo term (kondospectrumtk/twotime.py, edtwotimeref.py),
following Ternes, New J. Phys. 17 063016 (2015), arXiv:1505.04430.

This is a *different* way of computing the same third-order Kondo term
already covered (and validated against Fig. 3a/b) by
conductance.third_order_kondo_dIdV, which sums explicitly over
eigenstates. The two-time construction instead builds a Heisenberg
three-point function G(t2,tau)=<GS|Sl(t2+tau)Sk(t2)Sj(0)|GS> via real-time
evolution and extracts the same physical quantity through two closed-form
time-domain kernels (a Hilbert-transform-based Theta0 filter and a direct
K_W convolution) -- see twotime.py's module docstring for the full
derivation. It exists because, unlike the eigenstate-sum approach, it
generalizes to DMRG without ever diagonalizing/enumerating excited
states (the actual DMRG side reuses the same kernel machinery on top of
real TDVP time evolution instead of ED's eigenbasis-exact evolution).

The ED-based G(t2,tau) here (edtwotimeref.py) is deliberately built via
exact eigenbasis time evolution rather than any excited-state
diagonalization shortcut, so this test is a genuine end-to-end check of
the two-time/kernel pipeline itself, not just of ED's (already validated
elsewhere) ability to compute transition matrix elements directly.

These tests use deliberately modest grid resolution (fast: a few to ~20s
per test) and correspondingly loose tolerances -- finer grids converge
tighter (verified interactively during development: ~0.01% agreement at
higher resolution) but are too slow for routine test runs. See PR
history / kondospectrumtk module docstrings for the resolution/accuracy
tradeoffs (K_W needs t2 spacing finer than 1/omega0 and a range wider
than several/Gamma0; the Hilbert-transform-based Theta0 filter converges
much faster than that and is not the bottleneck).
"""
import numpy as np
import pytest

from dmrgpy import spinchain
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
from dmrgpy.kondospectrumtk.conductance import third_order_kondo_dIdV
from dmrgpy.kondospectrumtk.edtwotimeref import two_time_kondo_term_ed
from dmrgpy.kondospectrumtk.twotime import theta0_filter, K_W
from dmrgpy.kondospectrumtk.stepfunctions import Theta0, F0

G = 2.0
MUB = 5.7883818066e-5 # eV/T


def test_theta0_filter_reproduces_theta0_on_pure_exponentials():
    """Unit check of the Hilbert-transform-based Theta0 filter alone:
    for a pure exponential exp(-i*eps*tau) (a single "eigenstate" term),
    the filter should reproduce Theta0(eV-eps) essentially exactly (this
    is the piece that converges to machine precision even on coarse
    grids -- confirmed during development; the K_W piece is the one that
    needs a fine/wide grid)."""
    tau_grid = np.linspace(-2*np.pi/1e-5, 2*np.pi/1e-5, 4000, endpoint=False)
    eps = 3e-4
    eV = 5e-4
    signal = np.exp(-1j*eps*tau_grid)
    val = theta0_filter(tau_grid, signal, eV)
    assert val == pytest.approx(Theta0(np.array([eV-eps]))[0], abs=1e-6)


def test_K_W_matches_F0_pm_via_direct_integration():
    """Unit check of the K_W kernel alone: integrating K_W(t2;eV) against
    exp(-i*wm*t2) over t2 should reproduce F0(eV-wm)+F0(eV+wm)."""
    from scipy import integrate
    omega0, Gamma0 = 0.2, 5e-6
    eV, wm = 3e-4, 1.6e-4
    expected = F0(np.array([eV-wm]), omega0, Gamma0)[0] + F0(np.array([eV+wm]), omega0, Gamma0)[0]
    T = 30/Gamma0
    re, _ = integrate.quad(lambda t2: K_W(np.array([t2]), eV, omega0, Gamma0)[0]*np.cos(wm*t2),
                            -T, T, limit=2000, points=[0.])
    assert re == pytest.approx(expected, rel=2e-3)


def test_two_time_kondo_term_matches_eigenstate_sum():
    """End-to-end check: the two-time-correlator construction (real-time
    evolution + kernel extraction) must agree with the already-validated
    excited-state-sum third_order_kondo_dIdV for the same system, at
    T=0. Uses a small S=1/2 system and modest grid resolution (see module
    docstring) -- loose (but not vacuous) tolerance."""
    sc = spinchain.Spin_Chain(["1/2"])
    sc.set_hamiltonian(G*MUB*10.0*sc.Sz[0])
    ks = KondoSpectrum(sc, site=0, T=0.0)

    Jrho_s, T0 = -0.05, 1.0
    omega0, Gamma0 = 2e-3, 5e-6
    eVs = np.linspace(-2e-3, 2e-3, 9)

    mine = 4*np.pi*T0**2*Jrho_s*two_time_kondo_term_ed(
        ks, eVs, omega0=omega0, Gamma0=Gamma0,
        t2_width=25/Gamma0, t2_npts=40_000, t2_batch=10_000,
        tau_width=2*np.pi/2e-5, tau_npts=1_000)
    ref = third_order_kondo_dIdV(ks, eVs, Jrho_s, T0=T0, omega0=omega0, Gamma0=Gamma0)

    assert np.max(np.abs(mine-ref)) < 0.02
    assert np.max(np.abs(mine-ref)) < 0.05*np.max(np.abs(ref)) # also a relative check
