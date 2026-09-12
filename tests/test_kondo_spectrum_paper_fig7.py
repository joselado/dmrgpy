"""Pins of the STM/Kondo spectrum (kondospectrumtk/) against values
digitized from Fig. 7 of the arXiv v1 PDF of Ternes, New J. Phys. 17
063016 (2015), arXiv:1505.04430 (the paper's third-order S=1/2 figure:
panels b/d are absolutely scaled, in e^2 T0^2/h), plus regressions for
the three 2026-09-12 changes that those pins caught or motivated:

  - the potential-interference term's exchange diagram enters with the
    OPPOSITE sign to the direct one (F(eV-eps_im) - F(eV+eps_im)), so
    the term vanishes identically at B=0 and has no m=i zero-bias spike;
    the summed form it used to have put the 0 T, U=0.25 peak at -0.2 mV
    (Fig. 7d has it at 0) and made the 10 T step asymmetry ~0.27
    against the figure's ~0.05;
  - F is the paper's own symmetric closed form (eq. 22), which puts the
    Fig. 7b tails at 0.885 -- the electron-like eq. 20 used for both
    diagrams before gave 0.854;
  - the third-order sums run over thermally occupied initial states
    only, and T=0 averages over a degenerate ground-state manifold.

The digitized values come from the curves' own pixel colours on a
600 dpi render, with the axes calibrated on the tick marks (the 7d panel
runs to 1.45, not 1.4); the dashed second-order references in the same
panels then read 0.750/0.250 (7b) and 1.000/0.500 (7d), so the pins are
good to ~0.005. Comparing at the 10 T tails and step overshoots tests
the odd term's sign AND magnitude, which the zero-bias values (where it
vanishes) cannot.
"""
import numpy as np
import pytest

from dmrgpy import spinchain
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
from dmrgpy.kondospectrumtk.conductance import (
    third_order_kondo_dIdV, third_order_potential_dIdV, _occupied_states)
from dmrgpy.kondospectrumtk.stepfunctions import FBuilder, Theta, F0

G = 2.0
MUB = 5.7883818066e-5 # eV/T
KB = 8.617333262e-5 # eV/K
TP = 2*np.pi # get_kondo_spectrum returns 2*pi times the figure's units


def _chain(B):
    sc = spinchain.Spin_Chain(["1/2"])
    sc.set_hamiltonian(G*MUB*B*sc.Sz[0])
    return sc


def _spectrum(B, U, eVs):
    _, d = _chain(B).get_kondo_spectrum(eVs, site=0, Jrho_s=-0.05, U=U,
                                        T=1.0, order=3, omega0=20e-3)
    return d/TP


def test_fig7b_zero_bias_and_tails():
    """U=0: zero-bias values at 0, 0.5, 1, 2.5, 10 T and the +-4 mV
    tails, which every field curve shares."""
    eVs = np.array([-4e-3, 0., 4e-3])
    for B, y0 in ((0., 1.13), (0.5, 1.082), (1., 0.942), (2.5, 0.535),
                  (10., 0.324)):
        d = _spectrum(B, 0.0, eVs)
        assert d[1] == pytest.approx(y0, abs=0.01)
        assert d[0] == pytest.approx(0.886, abs=0.006) # tails
        assert d[2] == pytest.approx(0.886, abs=0.006)


def test_fig7d_zero_bias_values_and_symmetric_zero_field_peak():
    """U=0.25: the odd term vanishes at eV=0, so every zero-bias value is
    the U=0 one plus the elastic 4U^2=0.25 -- 1.39, 1.332, 1.192, 0.786,
    0.574 -- and at B=0 the whole curve is exactly symmetric with its
    peak at eV=0 (the summed exchange sign put it at -0.2 mV, 1.47)."""
    eVs = np.linspace(-4e-3, 4e-3, 161)
    for B, y0 in ((0., 1.39), (0.5, 1.332), (1., 1.192), (2.5, 0.786),
                  (10., 0.574)):
        d = _spectrum(B, 0.25, eVs)
        assert d[80] == pytest.approx(y0, abs=0.01)
    d0 = _spectrum(0., 0.25, eVs)
    assert np.allclose(d0, d0[::-1], atol=1e-10)
    assert np.argmax(d0) == 80
    assert d0[0] == pytest.approx(1.135, abs=0.006)


def test_fig7d_10T_asymmetry_matches_the_figure():
    """The odd (potential-interference) term's sign and size: at 10 T the
    figure's step overshoots read 1.231 (-1.67 mV) / 1.177 (+1.69 mV) and
    the +-4 mV tails 1.146 / 1.128."""
    eVs = np.linspace(-4e-3, 4e-3, 801)
    d = _spectrum(10., 0.25, eVs)
    left, right = (eVs < -1e-3), (eVs > 1e-3)
    assert d[left].max() == pytest.approx(1.231, abs=0.012)
    assert d[right].max() == pytest.approx(1.177, abs=0.012)
    assert d[0] == pytest.approx(1.146, abs=0.008)
    assert d[-1] == pytest.approx(1.128, abs=0.008)
    assert d[left].max() - d[right].max() == pytest.approx(0.054, abs=0.012)
    assert d[0] - d[-1] == pytest.approx(0.018, abs=0.006)


def test_potential_term_vanishes_at_zero_field_at_any_T():
    """Direct minus exchange: with every eps_im=0 the bracket is 0."""
    eVs = np.linspace(-3e-3, 3e-3, 31)
    for T in (0., 1.0):
        ks = KondoSpectrum(_chain(0.), site=0, T=T)
        dU = third_order_potential_dIdV(ks, eVs, -0.05, 0.25, T0=1.0)
        assert np.allclose(dU, 0., atol=1e-14)
        # while the Kondo term is a real peak there
        assert third_order_kondo_dIdV(ks, np.array([0.]), -0.05)[0] > 1.


def test_potential_term_exchange_diagram_sign_against_hand_built_sum():
    """The (121u)/(121uR) structure of the paper's Fig. 7c, written out
    by hand for the S=1/2 doublet: c*[g(eV)-g(-eV)] with
    g(v) = Theta(v) * 1/2 * [F(v-D) - F(v+D)] -- the m=i loop (weight
    1/4, eps=0) drops out, and the reversed order carries a minus sign."""
    B = 10.
    D = G*MUB*B
    T = 1.0; kT = KB*T
    ks = KondoSpectrum(_chain(B), site=0, T=T)
    Fb = FBuilder(T)
    eVs = np.linspace(-3e-3, 3e-3, 25)
    got = third_order_potential_dIdV(ks, eVs, -0.05, 0.25, T0=1.0, Fb=Fb)
    def g(v):
        return Theta(v/kT)*0.5*(Fb(v - D) - Fb(v + D))
    ref = 4*np.pi*(-0.05)*0.25*(g(eVs) - g(-eVs))
    assert np.allclose(got, ref, atol=1e-8)


def test_occupied_state_restriction_is_exact():
    """Restricting the initial-state sum to p_i > P_CUT*max(p) changes
    nothing measurable: compare against p_cut=0 (every state) at a T where
    several states are occupied, on a system with a real spectrum."""
    sc = spinchain.Spin_Chain(["1", "1/2"])
    h = 3e-4*sc.Sz[0]*sc.Sz[0] + G*MUB*2.0*sc.Sz[0] + 2e-4*(
        sc.Sx[0]*sc.Sx[1] + sc.Sy[0]*sc.Sy[1] + sc.Sz[0]*sc.Sz[1])
    sc.set_hamiltonian(h)
    T = 0.3 # kT = 26 ueV: the ~1 meV spread of the six levels then spans
            # ~40 kT, so the top ones fall below P_CUT while several stay
    ks = KondoSpectrum(sc, site=0, T=T)
    assert 1 < len(_occupied_states(ks)) < ks.dim
    Fb = FBuilder(T)
    eVs = np.linspace(-2e-3, 2e-3, 9)
    for f in (third_order_kondo_dIdV,):
        a = f(ks, eVs, -0.05, Fb=Fb)
        b = f(ks, eVs, -0.05, Fb=Fb, p_cut=0.)
        assert np.allclose(a, b, rtol=1e-9, atol=1e-12)
    a = third_order_potential_dIdV(ks, eVs, -0.05, 0.2, Fb=Fb)
    b = third_order_potential_dIdV(ks, eVs, -0.05, 0.2, Fb=Fb, p_cut=0.)
    assert np.allclose(a, b, rtol=1e-9, atol=1e-12)
    # and at T=0 exactly one state is kept
    assert len(_occupied_states(KondoSpectrum(sc, site=0, T=0.))) == 1


def test_degenerate_ground_state_at_T0_is_averaged():
    """S=1 with D*Sz^2 + g*muB*B*Sz at the |0>/|-1> crossing g*muB*B=D:
    two accidentally degenerate ground states with different excitation
    spectra. T=0 must give the T->0+ limit (equal weights), which is also
    what a tiny T gives -- not whichever basis eigh returned."""
    D = 5e-4
    sc = spinchain.Spin_Chain(["1"])
    sc.set_hamiltonian(D*sc.Sz[0]*sc.Sz[0] + D*sc.Sz[0])
    ks0 = KondoSpectrum(sc, site=0, T=0.)
    assert np.allclose(ks0.p[:2], 0.5) and ks0.p[2] == 0.
    eVs = np.linspace(-1.5e-3, 1.5e-3, 13)
    d0 = third_order_kondo_dIdV(ks0, eVs, -0.05)
    # the two members of the manifold individually give different spectra
    ks_a = KondoSpectrum(sc, site=0, T=0.); ks_a.p = np.array([1., 0., 0.])
    ks_b = KondoSpectrum(sc, site=0, T=0.); ks_b.p = np.array([0., 1., 0.])
    da, db = (third_order_kondo_dIdV(k, eVs, -0.05) for k in (ks_a, ks_b))
    assert not np.allclose(da, db, atol=1e-6)
    assert np.allclose(d0, 0.5*(da + db), atol=1e-12)
    # and a tiny T (4e-3 K: kT = 3.5e-7 eV, far below D but the F/Theta
    # broadening is then also tiny) agrees with T=0 to the broadening
    ksT = KondoSpectrum(sc, site=0, T=4e-3)
    dT = third_order_kondo_dIdV(ksT, eVs, -0.05)
    assert np.allclose(dT, d0, atol=0.02*np.max(np.abs(d0)))


def test_FBuilder_table_matches_direct_convolution_and_T0_limit():
    # the residual is the convolution's own Simpson error at the cusp the
    # singularity subtraction leaves (h^2-scaling: 2e-8 at 1 K, 2e-5 at
    # 20 K where the u-grid step is 20x larger), not the spline's
    for T, w0 in ((1.0, 20e-3), (20.0, 20e-3), (0.1, 0.2)):
        fb = FBuilder(T, omega0=w0)
        xs = np.concatenate([np.linspace(-3*w0, 3*w0, 301),
                             np.linspace(-30*fb.kT, 30*fb.kT, 301)])
        assert np.allclose(fb(xs), fb._direct(xs), rtol=0., atol=5e-5)
        assert np.allclose(fb(xs), fb(-xs), rtol=0., atol=1e-9) # even
    # T -> 0 recovers the closed form F0 away from its Gamma0 peak
    fb = FBuilder(1e-3, omega0=20e-3)
    xs = np.array([-5e-3, -5e-4, 5e-4, 5e-3])
    assert np.allclose(fb(xs), F0(xs, 20e-3), rtol=1e-4)
