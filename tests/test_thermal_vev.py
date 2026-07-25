"""Regression coverage for Many_Body_Chain.vev(mode="ED", T=...) (finite-
temperature expectation values via edtk/edchain.py -> vevtk/thermalvev.py).

thermal_vev_ex() sums Boltzmann-weighted <n|O|n> over a set of exact
eigenstates obtained from algebra.lowest_states(), whose own default is
n=10 -- previously that default was used unconditionally, so any chain
with a Hilbert space bigger than 10 states silently summed over an
arbitrary, temperature-independent 10-state subset instead of the full
spectrum, giving a wrong thermal average with no error raised (only a
weak, easy-to-miss "highest excited state comparable to thermal energy"
print). Fixed by (1) defaulting n to the full Hilbert-space dimension
whenever that's small enough for exact diagonalization anyway (matching
what algebra.lowest_states does internally regardless), and (2) raising
instead of merely printing when the states actually used still carry
non-negligible Boltzmann weight -- see edtk/edchain.py::EDchain.vev and
vevtk/thermalvev.py::thermal_vev_ex.
"""
import numpy as np
import pytest

from dmrgpy import spinchain


def test_single_spin_thermal_magnetization_matches_analytic():
    """A single S=1/2 spin in a field, H = B*Sz, has exact eigenvalues
    +-B/2 with <Sz> = +-1/2, so the thermal average has the textbook
    closed form <Sz> = -1/2 * tanh(beta*B/2) -- an independent analytic
    check, not a regression-pinned value."""
    sc = spinchain.Spin_Chain(["S=1/2"])
    B = 1.3
    sc.set_hamiltonian(B * sc.Sz[0])

    for T in (0.3, 0.7, 2.0, 10.0):
        beta = 1. / T
        expected = -0.5 * np.tanh(beta * B / 2.)
        got = sc.vev(sc.Sz[0], mode="ED", T=T)
        assert got.real == pytest.approx(expected, abs=1e-8)
        assert got.imag == pytest.approx(0., abs=1e-8)


def test_thermal_vev_default_matches_full_spectrum_reference():
    """6-site Heisenberg chain + field (Hilbert dim 64): before the fix,
    the default call silently summed over only the 10 lowest eigenstates
    (algebra.lowest_states' own default) regardless of dim or T, which at
    T=2.0 here previously returned ~-0.0975 instead of the true thermal
    average. The chain's Hilbert space is small enough to diagonalize
    exactly, so the default must now match an independently computed,
    full-spectrum Boltzmann sum built directly with numpy.linalg.eigh
    (bypassing thermal_vev_ex entirely) to any numerical precision."""
    n = 6
    sc = spinchain.Spin_Chain(["S=1/2" for _ in range(n)])
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]
    for i in range(n):
        h = h + 0.3 * sc.Sz[i]
    sc.set_hamiltonian(h)

    T = 2.0
    edobj = sc.get_ED_obj()
    Hd = np.array(edobj.get_hamiltonian().todense())
    assert Hd.shape[0] == 2 ** n

    es, vs = np.linalg.eigh(Hd)
    es = es - es.min()
    P = np.exp(-es / T)
    P /= P.sum()
    Opd = np.array(edobj.MO2matrix(sc.Sz[0]).todense())
    vals = np.array([np.conjugate(vs[:, i]) @ Opd @ vs[:, i] for i in range(len(es))])
    expected = np.sum(P * vals).real

    got = sc.vev(sc.Sz[0], mode="ED", T=T)
    assert got.real == pytest.approx(expected, abs=1e-8)

    # explicitly requesting the full spectrum must agree too
    got_explicit = sc.vev(sc.Sz[0], mode="ED", T=T, n=Hd.shape[0])
    assert got_explicit.real == pytest.approx(expected, abs=1e-8)


def test_thermal_vev_raises_instead_of_silently_truncating():
    """If the caller explicitly requests too few excited states for the
    requested temperature, thermal_vev_ex must raise rather than
    silently return a wrong average (the previous behavior: a print
    warning, then a wrong number)."""
    n = 6
    sc = spinchain.Spin_Chain(["S=1/2" for _ in range(n)])
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]
    for i in range(n):
        h = h + 0.3 * sc.Sz[i]
    sc.set_hamiltonian(h)

    with pytest.raises(RuntimeError):
        sc.vev(sc.Sz[0], mode="ED", T=2.0, n=5)

    # low enough T that only a handful of states carry any weight: same
    # small n must NOT raise, and must still match the full-spectrum value
    edobj = sc.get_ED_obj()
    Hd = np.array(edobj.get_hamiltonian().todense())
    es, vs = np.linalg.eigh(Hd)
    es = es - es.min()
    T_low = 0.05
    P = np.exp(-es / T_low)
    P /= P.sum()
    Opd = np.array(edobj.MO2matrix(sc.Sz[0]).todense())
    vals = np.array([np.conjugate(vs[:, i]) @ Opd @ vs[:, i] for i in range(len(es))])
    expected_low = np.sum(P * vals).real

    got_low = sc.vev(sc.Sz[0], mode="ED", T=T_low, n=5)
    assert got_low.real == pytest.approx(expected_low, abs=1e-6)
