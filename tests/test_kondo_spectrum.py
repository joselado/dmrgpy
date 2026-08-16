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


def test_third_order_kondo_zero_field_is_a_symmetric_zero_bias_peak():
    """Fig. 3a/b: at zero field the third-order Kondo term is a resonance
    centered exactly at eV=0 and symmetric in bias. Both tunneling
    directions contribute with the same sign (eq. "sym_z"), which is what
    makes it even -- every per-process curve of Fig. 3a is drawn symmetric
    about zero bias, which d(I^{t->s})/dV alone (a one-sided Theta step)
    cannot be. Regression test for a real bug: the s->t direction used to
    be omitted, leaving a one-sided half-peak whose maximum sat at the
    first positive bias point instead of at eV=0."""
    ks = _single_spin_half(B=0.0, T=1.0)
    span = 3e-3
    eVs = np.linspace(-span, span, 41) # odd count -> eVs[20] is exactly 0
    d3 = third_order_kondo_dIdV(ks, eVs, Jrho_s=-0.05, T0=1.0)
    assert np.all(third_order_kondo_dIdV(ks, eVs, Jrho_s=0., T0=1.0) == 0.)
    assert np.allclose(d3, d3[::-1], atol=1e-12) # even in eV
    assert np.argmax(d3) == len(eVs)//2 # peak exactly at eV=0
    assert d3[len(eVs)//2] > 2*d3[0] # a real peak, not a flat background


def test_third_order_potential_vanishes_at_U0_and_is_odd_in_bias():
    """Fig. 3c: the potential-interference term's two tunneling directions
    enter with opposite signs (eq. "asym_U"), so the term is odd in eV --
    the "sum" curve of Fig. 3c is antisymmetric about zero bias and
    crosses through it. Regression test for the same omitted-s->t bug as
    above: this term used to be one-sided (identically zero at negative
    bias for T=0, and finite at eV=0), which is not merely asymmetric but
    the wrong shape entirely."""
    ks = _single_spin_half(B=10.0, T=1.0)
    eVs = np.linspace(-2e-3, 2e-3, 21) # odd count -> eVs[10] is exactly 0
    dU0 = third_order_potential_dIdV(ks, eVs, Jrho_s=-0.05, U=0.0, T0=1.0)
    assert np.allclose(dU0, 0.)
    dU = third_order_potential_dIdV(ks, eVs, Jrho_s=-0.05, U=0.25, T0=1.0)
    assert not np.allclose(dU, dU[::-1], atol=1e-8) # bias asymmetric
    assert np.allclose(dU, -dU[::-1], atol=1e-12) # and specifically odd
    assert dU[len(eVs)//2] == pytest.approx(0., abs=1e-14)
    assert np.max(np.abs(dU)) > 1e-3 # not trivially zero everywhere


def test_third_order_matches_paper_figure_3b_and_3d_peak_values():
    """Absolute cross-check of every prefactor in the third-order spectrum
    against the paper's own Figs. 3b/3d, which are plotted in absolute
    units (e^2 T0^2/h) for a single S=1/2 at B=0, T=1K, Jrho_s=-0.05 and
    the paper-default omega0=20meV. This function's own return convention
    is 2*pi times those figure units (its second-order plateau is 2*pi*0.75
    where Fig. 3b's dashed curve reads 0.75).

    Reading the zero-bias peaks off the published figures gives 1.13
    (Fig. 3b, U=0) and 1.39 (Fig. 3d, U=0.25). These pin down three
    separate normalizations that no other test here constrains: the
    Levi-Civita coefficient Im[X]/2 (an Im[X]/4 reading of eq.
    "3rd-normal" misses a factor 2), the sum over both tunneling
    directions, and the elastic potential channel's 4*|U|^2 (eq.
    "Matrix1sq"'s bare |U|^2 misses a factor 4 -- Fig. 3d minus Fig. 3b
    for the dashed second-order curves is 0.25 = 4*0.25^2). The
    potential-interference term contributes exactly 0 at eV=0 (it is odd),
    so 3d - 3b at zero bias is the elastic term alone."""
    ks = _single_spin_half(B=0.0, T=1.0)
    omega0 = 20e-3 # "If not otherwise noted we use omega0=20 meV"
    z = np.array([0.0])
    Jrho_s, U, tp = -0.05, 0.25, 2*np.pi

    kondo = third_order_kondo_dIdV(ks, z, Jrho_s, T0=1.0, omega0=omega0)[0]
    pot = third_order_potential_dIdV(ks, z, Jrho_s, U, T0=1.0, omega0=omega0)[0]
    fig3b = (second_order_dIdV(ks, z, T0=1.0, U=0.0)[0] + kondo)/tp
    fig3d = (second_order_dIdV(ks, z, T0=1.0, U=U)[0] + kondo + pot)/tp

    assert fig3b == pytest.approx(1.13, abs=0.02)
    assert fig3d == pytest.approx(1.39, abs=0.02)
    # the elastic potential channel alone, Fig. 3d - Fig. 3b (dashed)
    elastic = (second_order_dIdV(ks, z, T0=1.0, U=U)[0]
               - second_order_dIdV(ks, z, T0=1.0, U=0.0)[0])/tp
    assert elastic == pytest.approx(4*U**2, abs=1e-12)
    assert elastic == pytest.approx(0.25, abs=1e-12)


def test_get_kondo_spectrum_order2_matches_second_order_dIdV():
    sc = spinchain.Spin_Chain(["1/2"])
    sc.set_hamiltonian(G*MUB*5.0*sc.Sz[0])
    eVs = np.linspace(-1e-3, 1e-3, 11)
    eV_out, dIdV = sc.get_kondo_spectrum(eVs, site=0, T=1.0, order=2)
    ks = KondoSpectrum(sc, site=0, T=1.0)
    ref = second_order_dIdV(ks, eVs, T0=1.0, U=0.0)
    assert np.allclose(eV_out, eVs)
    assert np.allclose(dIdV, ref)


def test_third_order_kondo_eps_if_sign_on_asymmetric_spectrum():
    """Regression test for a real sign bug: third_order_kondo_dIdV used
    to build eps_if as e_i-e_f instead of e_f-e_i (the convention its own
    module docstring and second_order_dIdV both state), silently -- every
    other test/example here only exercises symmetric two-level (S=1/2)
    systems, where swapping i<->f just permutes which pair gets which
    energy and the bug has no effect on the Boltzmann-weighted sum. An
    S=1 impurity with both anisotropy and a field has a non-degenerate,
    asymmetric-under-relabeling spectrum, so it does distinguish the two
    conventions: cross-check third_order_kondo_dIdV against an
    independently written reference that only uses the (m,i) energy
    differences and Boltzmann weights, not the module's own eps_if."""
    sc = spinchain.Spin_Chain(["1"])
    D, B = 3e-4, 4.0
    sc.set_hamiltonian(D*sc.Sz[0]*sc.Sz[0] + G*MUB*B*sc.Sz[0])
    ks = KondoSpectrum(sc, site=0, T=1.0)
    assert len(set(np.round(ks.e, 12))) == ks.dim # genuinely non-degenerate

    eVs = np.linspace(-2e-3, 2e-3, 11)
    Jrho_s = -0.05
    got = third_order_kondo_dIdV(ks, eVs, Jrho_s, T0=1.0)

    # independent reference: same physics, written from scratch here
    Xi = np.stack([ks.Sx, ks.Sy, ks.Sz], axis=-1)
    eps3 = np.zeros((3, 3, 3))
    for a, b, c in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
        eps3[a, b, c] = 1.
    for a, b, c in [(0, 2, 1), (2, 1, 0), (1, 0, 2)]:
        eps3[a, b, c] = -1.
    coeff = np.imag(np.einsum('jkl,ifl,fmk,mij->ifm', eps3, Xi, Xi, Xi))/2.
    kT = ks.kB*ks.T
    from dmrgpy.kondospectrumtk.stepfunctions import FBuilder
    Fb = FBuilder(ks.T) # build once, reuse (F() alone would be far too slow here)
    ref = np.zeros(len(eVs))
    for e_idx, eV in enumerate(eVs):
        acc = 0.
        # both tunneling directions: v=+eV is t->s, v=-eV is s->t, added
        # with the same sign for the Kondo term (eq. "sym_z")
        for v in (eV, -eV):
            for i in range(ks.dim):
                for f in range(ks.dim):
                    eps_if_ref = ks.e[f] - ks.e[i] # unambiguous by construction
                    th = Theta(np.array([(v - eps_if_ref)/kT]))[0]
                    for m in range(ks.dim):
                        eps_im_ref = ks.e[m] - ks.e[i]
                        fsum = (Fb(np.array([v - eps_im_ref]))[0]
                                + Fb(np.array([v + eps_im_ref]))[0])
                        acc += ks.p[i]*coeff[i, f, m]*th*fsum
        ref[e_idx] = 4*np.pi*1.0**2*Jrho_s*acc
    assert np.allclose(got, ref, atol=1e-8)


def test_get_kondo_spectrum_order3_is_finite_and_real():
    sc = spinchain.Spin_Chain(["1/2"])
    sc.set_hamiltonian(G*MUB*5.0*sc.Sz[0])
    eVs = np.linspace(-2e-3, 2e-3, 21)
    eV_out, dIdV = sc.get_kondo_spectrum(
            eVs, site=0, Jrho_s=-0.05, U=0.2, T=1.0, order=3)
    assert dIdV.dtype == np.float64
    assert np.all(np.isfinite(dIdV))


def test_get_kondo_spectrum_mode_ed_T0_matches_direct_call():
    """mode="ED" (the default) accepts T=0 directly, matching a manual
    KondoSpectrum(T=0) + conductance call."""
    sc = spinchain.Spin_Chain(["1/2"])
    sc.set_hamiltonian(G*MUB*5.0*sc.Sz[0])
    eVs = np.linspace(-1e-3, 1e-3, 11)
    eV_out, dIdV = sc.get_kondo_spectrum(eVs, site=0, Jrho_s=-0.05, T=0.0,
                                          order=3)
    ks = KondoSpectrum(sc, site=0, T=0.0)
    ref = (second_order_dIdV(ks, eVs, T0=1.0, U=0.0)
           + third_order_kondo_dIdV(ks, eVs, -0.05, T0=1.0))
    assert np.allclose(dIdV, ref)


def test_get_kondo_spectrum_mode_dmrg_rejects_nonzero_T_and_missing_kwargs():
    sc = spinchain.Spin_Chain(["1/2"])
    sc.set_hamiltonian(G*MUB*5.0*sc.Sz[0])
    eVs = np.linspace(-1e-3, 1e-3, 5)
    with pytest.raises(ValueError):
        sc.get_kondo_spectrum(eVs, site=0, T=1.0, mode="DMRG")
    # U!=0, order=3 is implemented for mode="DMRG" (potentialdc.py), but
    # still requires the same mandatory `es` kwarg as the plain
    # second-order term (checked before the potential term is reached) --
    # see test_kondo_spectrum_potentialdc.py and
    # test_kondo_spectrum_dmrgtwotime.py for the functional coverage of
    # this path with a compiled backend.
    with pytest.raises(ValueError):
        sc.get_kondo_spectrum(eVs, site=0, T=0.0, U=0.1, order=3, mode="DMRG")
    with pytest.raises(ValueError):
        sc.get_kondo_spectrum(eVs, site=0, mode="bogus")


def test_get_kondo_spectrum_dmrg_combines_terms_correctly(monkeypatch):
    """Unit-level check of Spin_Chain._get_kondo_spectrum_dmrg's wiring --
    the correct linear combination of the second-order term, third-order
    Kondo term, and (if U!=0) the potential-interference term -- with the
    three underlying (expensive) computations monkeypatched out entirely.
    The individual pieces are numerically validated against their ED
    references elsewhere (test_kondo_spectrum_potentialdc.py,
    test_kondo_spectrum_dmrgtwotime.py); a real end-to-end run through
    mode="DMRG" submode="KPM" was tried here first and found
    impractically expensive for a routine test -- a single KPM dynamical-
    correlator call took ~40s on this repo's test hardware even with a
    trivially small number of moments, confirmed directly -- so this
    checks the arithmetic only, with no real DMRG computation and no
    compiled-backend dependency (the local `from .kondospectrumtk.X
    import Y` imports inside _get_kondo_spectrum_dmrg re-resolve Y from
    X's namespace on every call, so patching X.Y here is picked up)."""
    import dmrgpy.kondospectrumtk.secondorder_dc as secondorder_dc
    import dmrgpy.kondospectrumtk.potentialdc as potentialdc
    import dmrgpy.kondospectrumtk.dmrgtwotime as dmrgtwotime

    eVs = np.linspace(-1e-3, 1e-3, 3)
    monkeypatch.setattr(secondorder_dc, "second_order_dIdV_dc",
                         lambda *a, **kw: np.full(len(eVs), 1.0))
    monkeypatch.setattr(dmrgtwotime, "two_time_kondo_term_dmrg",
                         lambda *a, **kw: np.full(len(eVs), 2.0))
    monkeypatch.setattr(potentialdc, "third_order_potential_dIdV_dc",
                         lambda *a, **kw: np.full(len(eVs), 3.0))

    sc = spinchain.Spin_Chain(["1/2"])
    sc.set_hamiltonian(G*MUB*5.0*sc.Sz[0])
    T0, Jrho_s = 2.0, 0.5
    common = dict(es=np.array([0.]), dt2=1., n_t2_half=1, dtau=1.,
                  n_tau_half=1)

    _, dIdV = sc._get_kondo_spectrum_dmrg(eVs, 0, Jrho_s, 0.0, T0, 20e-3,
                                           5e-6, order=3, **common)
    assert np.allclose(dIdV, 1.0 + 4*np.pi*T0**2*Jrho_s*2.0) # U=0: no potential term

    _, dIdV = sc._get_kondo_spectrum_dmrg(eVs, 0, Jrho_s, 0.3, T0, 20e-3,
                                           5e-6, order=3, **common)
    assert np.allclose(dIdV, 1.0 + 4*np.pi*T0**2*Jrho_s*2.0 + 3.0)

    _, dIdV = sc._get_kondo_spectrum_dmrg(eVs, 0, Jrho_s, 0.3, T0, 20e-3,
                                           5e-6, order=2, **common)
    assert np.allclose(dIdV, 1.0) # order=2: no third-order terms at all


def test_two_time_kondo_term_dmrg_requires_explicit_grid():
    """dt2/n_t2_half/dtau/n_tau_half have no numeric default (see
    two_time_kondo_term_dmrg's docstring in dmrgtwotime.py): a small, fast
    default silently returns a badly wrong result for the also-default
    omega0/Gamma0 instead of erroring, so all four are mandatory. This
    check fires before the (here, nonsense) chain/site arguments are ever
    touched, so it needs no compiled itensor_version=3 backend -- unlike
    the rest of dmrgtwotime.py's coverage in
    test_kondo_spectrum_dmrgtwotime.py, which is skipped without one."""
    from dmrgpy.kondospectrumtk.dmrgtwotime import two_time_kondo_term_dmrg
    with pytest.raises(ValueError):
        two_time_kondo_term_dmrg(None, 0, [0.0])
    with pytest.raises(ValueError):
        two_time_kondo_term_dmrg(None, 0, [0.0], dt2=1.0, n_t2_half=5,
                                  dtau=1.0, n_tau_half=None)
