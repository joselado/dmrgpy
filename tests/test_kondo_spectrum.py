"""Coverage for the third-order STM/Kondo perturbation-theory spectrum
(kondospectrumtk/, Spin_Chain.get_kondo_spectrum), following Ternes, New
J. Phys. 17 063016 (2015), arXiv:1505.04430.

These checks are analytic/structural, not golden-value regressions: for
a single S=1/2 impurity the second-order term's plateau values can be
worked out exactly by hand (see second_order_dIdV's own derivation in
kondospectrumtk/conductance.py), and several qualitative features (the
paper's own Fig. 2/3) are unambiguous enough to assert on directly:
Theta(x) must be a bounded, symmetric step; F(eps,T) must be a positive,
peaked, decaying function (checked against digitized-by-eye values from
the paper's own Fig. F); the zero-field third-order Kondo term must be a
peak at eV=0; the potential-interference term must vanish at U=0 and be
bias-asymmetric at U!=0.
"""
import numpy as np
import pytest

from dmrgpy import spinchain
from dmrgpy.kondospectrumtk.stepfunctions import Theta, Theta_prime, F
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
from dmrgpy.kondospectrumtk.conductance import (
    second_order_dIdV, third_order_kondo_dIdV, third_order_potential_dIdV)

G = 2.0
MUB = 5.7883818066e-5 # eV/T
KB = 8.617333262e-5 # eV/K


def _single_spin_half(B, T=1.0):
    sc = spinchain.Spin_Chain(["1/2"])
    sc.set_hamiltonian(G*MUB*B*sc.Sz[0])
    return KondoSpectrum(sc, site=0, T=T)


def test_theta_bounded_and_symmetric():
    xs = np.linspace(-30, 30, 61)
    th = Theta(xs)
    assert np.all(th >= 0.) and np.all(th <= 1.)
    assert np.allclose(th + Theta(-xs), 1.)
    assert Theta(np.array([0.]))[0] == pytest.approx(0.5)


def test_theta_prime_matches_finite_difference():
    xs = np.linspace(-5, 5, 21)
    h = 1e-6
    fd = (Theta(xs+h) - Theta(xs-h))/(2*h)
    assert np.allclose(Theta_prime(xs), fd, atol=1e-4)


def test_F_matches_figure_F_reference_values():
    # digitized from the paper's own Fig. F(b), T=1K, omega0=200meV,
    # Gamma0=5ueV curve: F(0)~7.5, F(10meV)~3.1
    vals = F(np.array([0., 10e-3]), T=1.0, omega0=0.2, Gamma0=5e-6)
    assert vals[0] == pytest.approx(7.5, abs=0.1)
    assert vals[1] == pytest.approx(3.1, abs=0.2)
    assert vals[0] > vals[1] > 0. # peaked and positive


def test_F_decays_away_from_the_peak():
    # F is "electron-like" (equ. "F_1" uses 1-f(ep',T) in its numerator,
    # not the symmetrized combination with the "hole-like" equ. "F_1h"),
    # so it is NOT expected to be even in x=eV-eps_m -- only decay away
    # from its peak at x=0 in both directions.
    xs = np.array([0., 1e-3, 1.])
    vals = F(xs, T=1.0, omega0=0.2, Gamma0=5e-6)
    assert vals[0] > vals[1] > vals[2] # monotonically decaying away from 0
    assert vals[0] > 0. and vals[1] > 0. # positive near the peak
    neg_vals = F(np.array([-1e-3, -1.]), T=1.0, omega0=0.2, Gamma0=5e-6)
    assert vals[0] > neg_vals[0] > neg_vals[1]


def test_second_order_zeeman_step_plateaus():
    """Analytic plateaus for a single S=1/2 impurity, unpolarized leads,
    U=0, field large compared to kT (so the ground state is essentially
    fully occupied): at large positive bias only the t->s direction is
    above threshold for both its elastic (weight 0.25) and spin-flip
    (weight 0.5) channels, while s->t is deep below its own threshold
    (Theta(-eV-eps_if)~0) -- giving 2*pi*0.75. Exactly at eV=0 only the
    elastic channel sits at its own threshold in *both* directions at
    once (Theta(0)=1/2 each), giving 2*pi*2*(0.25*0.5) = pi/2."""
    ks = _single_spin_half(B=10.0, T=1.0)
    eVs = np.array([0., 5.0*G*MUB*10.0])
    dIdV = second_order_dIdV(ks, eVs, T0=1.0, U=0.0)
    assert dIdV[0] == pytest.approx(np.pi/2, abs=1e-3)
    assert dIdV[1] == pytest.approx(2*np.pi*0.75, abs=1e-3)


def test_second_order_symmetric_and_positive():
    ks = _single_spin_half(B=3.0, T=1.0)
    eVs = np.linspace(-1e-3, 1e-3, 41)
    dIdV = second_order_dIdV(ks, eVs, T0=1.0, U=0.15)
    assert np.all(dIdV > 0.)
    assert np.allclose(dIdV, dIdV[::-1], atol=1e-8)


def test_third_order_kondo_zero_field_has_a_zero_bias_feature():
    """Fig. 3a/b: at zero field the third-order Kondo term produces a
    zero-bias resonance. This term is d(I^{t->s})/dV only (see
    third_order_kondo_dIdV's docstring), so unlike the full, bidirectional
    Fig. 3b curve it is not itself symmetric in eV -- Theta(eV-eps_if) is
    a one-sided step -- so the peak sits strictly inside the swept window
    rather than exactly at eV=0, but should be far above both the Jrho_s=0
    baseline (zero) and the window edges."""
    ks = _single_spin_half(B=0.0, T=1.0)
    span = 3e-3
    eVs = np.linspace(-span, span, 41)
    d3 = third_order_kondo_dIdV(ks, eVs, Jrho_s=-0.05, T0=1.0)
    assert np.all(third_order_kondo_dIdV(ks, eVs, Jrho_s=0., T0=1.0) == 0.)
    imax = np.argmax(d3)
    assert 0 < imax < len(eVs) - 1 # peak strictly inside the window
    assert d3[imax] > d3[0] and d3[imax] > d3[-1]


def test_third_order_potential_vanishes_at_U0_and_is_asymmetric():
    ks = _single_spin_half(B=10.0, T=1.0)
    eVs = np.linspace(-2e-3, 2e-3, 21)
    dU0 = third_order_potential_dIdV(ks, eVs, Jrho_s=-0.05, U=0.0, T0=1.0)
    assert np.allclose(dU0, 0.)
    dU = third_order_potential_dIdV(ks, eVs, Jrho_s=-0.05, U=0.25, T0=1.0)
    assert not np.allclose(dU, dU[::-1], atol=1e-8) # bias asymmetric


def test_get_kondo_spectrum_order2_matches_second_order_dIdV():
    sc = spinchain.Spin_Chain(["1/2"])
    sc.set_hamiltonian(G*MUB*5.0*sc.Sz[0])
    eVs = np.linspace(-1e-3, 1e-3, 11)
    eV_out, dIdV = sc.get_kondo_spectrum(eVs, site=0, T=1.0, order=2)
    ks = KondoSpectrum(sc, site=0, T=1.0)
    ref = second_order_dIdV(ks, eVs, T0=1.0, U=0.0)
    assert np.allclose(eV_out, eVs)
    assert np.allclose(dIdV, ref)


def test_get_kondo_spectrum_order3_is_finite_and_real():
    sc = spinchain.Spin_Chain(["1/2"])
    sc.set_hamiltonian(G*MUB*5.0*sc.Sz[0])
    eVs = np.linspace(-2e-3, 2e-3, 21)
    eV_out, dIdV = sc.get_kondo_spectrum(
            eVs, site=0, Jrho_s=-0.05, U=0.2, T=1.0, order=3)
    assert dIdV.dtype == np.float64
    assert np.all(np.isfinite(dIdV))
