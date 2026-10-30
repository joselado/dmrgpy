"""Open and renormalized Wilson chains, Bruognolo et al., PRB 95, 121115(R)
(2017), arXiv:1611.05291 (wilsonchain.py).

The quadratic checks are exact single-particle statements: the open chain
reproduces the bath correlator of its star at every frequency, and the C2
renormalized chain reproduces the static correlator of the continuum at
every length, where the standard Wilson chain shows the mass-flow error.
The many-body checks put the dissipative oscillator and the spin-boson
model on a short C2 chain and compare DMRG on "python" with ED."""

import warnings

import numpy as np
import pytest

from dmrgpy import wilsonchain as wc

S = 0.4
SUBOHMIC = lambda e: e**S if e > 0.0 else 0.0
# asymmetric pseudogap, |e|^1/2 with weight 1 above and 1/2 below the Fermi
# level: Re G(0) = -int A/e = -(0.5/r)/(1.5/(r+1)) = -1 at r = 1/2
PSEUDOGAP = lambda e: abs(e)**0.5*(1.0 if e > 0.0 else 0.5)


@pytest.mark.parametrize("band,dos", [((0.0, 1.0), SUBOHMIC),
                                      ((-1.0, 1.0), PSEUDOGAP),
                                      ((-0.5, 1.0), None)])
def test_open_chain_is_exact(band, dos):
    # the chain with every fast bath and the last slow bath is an exact
    # rewriting of the star, at any frequency off the real axis
    n = 10
    owc = wc.open_wilson_chain(n, Lambda=2.0, band=band, dos=dos)
    x, w = wc.bath_star(Lambda=2.0, band=band, dos=dos, nint=n+15)
    for om in (0.3+0.01j, 1e-3+1e-5j, -0.2+0.1j, 2e-4+1e-6j):
        ref = np.sum(w/(om - x))
        assert abs(owc.bath_green(om)/ref - 1.0) < 1e-11


def test_star_keeps_the_static_correlator_of_the_continuum():
    x, w = wc.bath_star(Lambda=2.0, band=(0.0, 1.0), dos=SUBOHMIC, nint=30)
    assert np.sum(w/x) == pytest.approx((S+1.0)/S, rel=1e-12)


@pytest.mark.parametrize("band,dos,exact", [((0.0, 1.0), SUBOHMIC, -(S+1)/S),
                                            ((-1.0, 1.0), PSEUDOGAP, -1.0)])
def test_c2_chain_has_no_mass_flow(band, dos, exact):
    owc = wc.open_wilson_chain(20, Lambda=2.0, band=band, dos=dos)
    for n in (1, 2, 5, 10, 20):
        g = {k: wc.chain_green(*owc.renormalized(n, kind=k)).real
             for k in ("C0", "C1", "C2")}
        assert g["C2"] == pytest.approx(exact, abs=1e-10)
        assert abs(g["C0"] - exact) > 0.1 # no shifts: wrong at every length
    # the fast shifts alone get the static value back only for a long chain
    c1 = [wc.chain_green(*owc.renormalized(n, kind="C1")).real for n in (5, 20)]
    assert abs(c1[1] - exact) < 0.1*abs(c1[0] - exact)


def test_dissipative_oscillator_susceptibility():
    # chi = 1/(Omega - 2 alpha/s) for the continuum; the standard Wilson
    # chain drifts with its length and saturates at a wrong value
    alpha, Omega = 0.1, 1.0
    exact = 1.0/(Omega - 2.0*alpha/S)
    chi = {}
    for kind in ("standard", "C2"):
        chi[kind] = []
        for n in (4, 12, 40):
            timp, e, t = wc.boson_bath_chain(n, alpha=alpha, s=S, kind=kind)
            chi[kind].append(1.0/(Omega + timp**2*wc.chain_green(e, t).real))
    assert np.max(np.abs(np.array(chi["C2"]) - exact)) < 1e-10
    assert np.all(np.diff(chi["standard"]) > 0.0)
    assert abs(chi["standard"][-1] - exact) > 0.05


def test_symmetric_band_has_no_shifts():
    owc = wc.open_wilson_chain(20, Lambda=2.0, band=(-1.0, 1.0))
    assert np.max(np.abs(owc.eps)) < 1e-14
    assert np.max(np.abs(owc.dF)) < 1e-14
    assert np.max(np.abs(owc.dS)) < 1e-12
    # one factor Lambda per site
    assert np.allclose(owc.tS[5:]/owc.tS[4:-1], 0.5, rtol=1e-2)


def test_bath_spectra_carry_the_couplings():
    owc = wc.open_wilson_chain(6, Lambda=2.0, band=(0.0, 1.0), dos=SUBOHMIC)
    om = np.linspace(-1.0, 3.0, 40001)
    for k in range(6):
        g = owc.bath_spectrum(k, om, eta=1e-3)
        assert np.sum(g)*(om[1]-om[0]) == pytest.approx(owc.tF[k]**2, rel=1e-3)
    # the last slow bath lives below the last window, band*Lambda^-6
    z, r = owc.slow
    assert np.sum(r) == pytest.approx(1.0, abs=1e-14)
    assert np.all((z > 0.0) & (z <= 2.0**-6))


def test_bad_arguments_raise():
    with pytest.raises(ValueError):
        wc.open_wilson_chain(5, band=(0.1, 1.0))
    with pytest.raises(ValueError):
        wc.open_wilson_chain(5, nint=6)
    owc = wc.open_wilson_chain(5)
    with pytest.raises(ValueError):
        owc.renormalized(kind="C3")
    with pytest.raises(ValueError):
        owc.renormalized(n=6)


def test_oscillator_dmrg_against_ed_and_exact():
    # quadratic model on a C2 chain: <x> is linear in eps, and the
    # truncation to 8 bosons per site is the only error against the exact
    # single-particle chi
    bc = wc.oscillator_chain(3, alpha=0.1, s=S, Omega=1.0, eps=0.05, maxnb=8,
                             itensor_version="python")
    bc.maxm = 40
    bc.nsweeps = 20
    dm = -bc.vev(bc.x_imp).real/0.05
    ed = -bc.vev(bc.x_imp, mode="ED").real/0.05
    e, t = bc.chain
    exact = 1.0/(1.0 + bc.timp**2*wc.chain_green(e, t).real)
    assert dm == pytest.approx(ed, abs=1e-6)
    assert dm == pytest.approx(exact, abs=1e-5)
    assert exact == pytest.approx(2.0, abs=1e-10) # the continuum value


def test_spin_boson_dmrg_against_ed():
    bc = wc.spin_boson_chain(3, alpha=0.1, s=S, Delta=0.1, eps=0.02, maxnb=6,
                             itensor_version="python")
    bc.maxm = 40
    bc.nsweeps = 20
    assert bc.gs_energy() == pytest.approx(bc.gs_energy(mode="ED"), abs=1e-6)
    assert bc.vev(bc.x_imp).real == pytest.approx(
        bc.vev(bc.x_imp, mode="ED").real, abs=1e-6)


def test_anderson_chain_kinds():
    # a symmetric flat band gives the same chain on every renormalized kind,
    # and an asymmetric density of states shifts the onsite energies
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        e0 = [wc.anderson_chain(3, kind=k, mode="ED").gs_energy()
              for k in ("C0", "C2")]
        assert e0[0] == pytest.approx(e0[1], abs=1e-12)
        dos = lambda e: 1.0 + 0.8*e
        fc = wc.anderson_chain(3, kind="C2", dos=dos, itensor_version="python")
        fc.maxm = 40
        fc.nsweeps = 16
        assert fc.gs_energy() == pytest.approx(fc.gs_energy(mode="ED"), abs=1e-6)
        with pytest.raises(ValueError):
            wc.anderson_chain(3, kind="C5")
