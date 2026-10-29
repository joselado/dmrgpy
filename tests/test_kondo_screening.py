"""Wilson chain coefficients and the NRG-MPS screening diagnostics of
arXiv:2501.10746 (wilsonchain.py, kondoscreening.py)."""

import numpy as np
import pytest

from dmrgpy import cppext, wilsonchain


def _backend(version):
    marks = ()
    if version in (2, 3):
        marks = pytest.mark.skipif(not cppext.available(version),
            reason="needs the compiled itensor_version=%d extension" % version)
    return pytest.param(version, id=str(version), marks=marks)


@pytest.mark.parametrize("Lambda", [1.5, 2.0, 3.0, 10.0])
def test_flat_band_matches_wilson_closed_form(Lambda):
    eps, t = wilsonchain.wilson_chain(40, Lambda=Lambda)
    assert np.max(np.abs(eps)) < 1e-14
    assert np.max(np.abs(t/wilsonchain.flat_band_hoppings(40, Lambda)-1)) < 1e-12


def test_chain_is_a_unitary_image_of_the_star():
    xi, w = wilsonchain.star_levels(Lambda=2.0, band=(-1.0, 0.7), nint=8)
    eps, t = wilsonchain.wilson_chain(len(xi), Lambda=2.0, band=(-1.0, 0.7), nint=8)
    h = np.diag(eps) + np.diag(t, 1) + np.diag(t, -1)
    assert np.allclose(np.linalg.eigvalsh(h), np.sort(xi), atol=1e-13)


def test_a_dos_function_gives_the_same_chain_as_the_flat_default():
    eps, t = wilsonchain.wilson_chain(20, dos=lambda e: 3.0)
    assert np.max(np.abs(t/wilsonchain.flat_band_hoppings(20)-1)) < 1e-10


def _chain(version, B=0.0, n=5):
    fc = wilsonchain.anderson_chain(n, ed=-0.5, U=1.0, tD=0.3, Lambda=4.0,
                                    B=B, itensor_version=version)
    fc.maxm = 64
    fc.nsweeps = 16
    fc.cutoff = 1e-12
    return fc


@pytest.mark.parametrize("version", [_backend("python"), _backend(3)])
def test_screening_dmrg_against_ed(version):
    fc = _chain(version)
    ed = fc.get_kondo_screening(mode="ED")
    dm = fc.get_kondo_screening()
    assert dm["M_imp"] == pytest.approx(ed["M_imp"], abs=1e-6)
    assert np.max(np.abs(dm["M_n"]-ed["M_n"])) < 1e-6
    assert np.max(np.abs(dm["schmidt"]-ed["schmidt"])) < 1e-6


def test_singlet_identities():
    # M_env = -M_imp in any state with total Sz = 0, and the singlet makes
    # the single-occupancy block proportional to the identity
    r = _chain("python").get_kondo_screening(mode="ED")
    assert r["M_env"] == pytest.approx(-r["M_imp"], abs=1e-10)
    assert r["schmidt"] == pytest.approx([2**-0.5, 2**-0.5], abs=1e-10)
    assert 0.5 < r["M_imp"] < 1.0


def test_a_strong_field_polarizes_the_moment():
    r = _chain("python", B=0.2).get_kondo_screening(mode="ED")
    assert r["schmidt"][0] > 0.99
    # the bath feels the field too, and M_n is not a connected correlator,
    # so the bath's own polarization makes M_env positive
    assert r["M_env"] > 0.0


@pytest.mark.parametrize("version", [_backend("python"), _backend(3)])
def test_resonant_level_against_wick(version):
    # U=0: 4<Sz_0 Sz_j> = -2 |<c_0^dag c_j>|^2 and M_imp = 2 n (1-n)
    n, tD = 9, 0.3
    eps, t = wilsonchain.wilson_chain(n, Lambda=4.0)
    h = np.diag(np.concatenate([[0.0], eps]))
    h[0, 1] = h[1, 0] = tD
    for i in range(n-1):
        h[i+1, i+2] = h[i+2, i+1] = t[i]
    ev, U = np.linalg.eigh(h)
    G = U[:, ev < 0] @ U[:, ev < 0].T
    fc = wilsonchain.anderson_chain(n, ed=0.0, U=0.0, tD=tD, Lambda=4.0,
                                    itensor_version=version)
    fc.maxm = 100
    fc.nsweeps = 16
    fc.cutoff = 1e-12
    r = fc.get_kondo_screening()
    assert r["M_imp"] == pytest.approx(2*G[0, 0]*(1-G[0, 0]), abs=1e-5)
    assert np.max(np.abs(r["M_n"] + 2*G[0, 1:]**2)) < 1e-5


@pytest.mark.parametrize("version", [_backend("python"), _backend(3)])
def test_two_level_impurity_dmrg_against_ed(version):
    fc = wilsonchain.anderson_chain(3, ed=[-0.5, -0.5], U=1.0, tD=[0.3, 0.15],
                                    Lambda=4.0, itensor_version=version)
    fc.maxm = 64
    fc.nsweeps = 16
    fc.cutoff = 1e-12
    ed = fc.get_kondo_screening(mode="ED")
    dm = fc.get_kondo_screening()
    assert np.max(np.abs(np.array(dm["M_imp"])-ed["M_imp"])) < 1e-6
    assert np.max(np.abs(np.array(dm["M_env"])-np.array(ed["M_env"]))) < 1e-6
    assert np.max(np.abs(dm["schmidt"]-ed["schmidt"])) < 1e-6


def test_equal_couplings_give_an_su4_singlet():
    fc = wilsonchain.anderson_chain(3, ed=[-0.5, -0.5], U=1.0, tD=[0.3, 0.3],
                                    Lambda=4.0)
    r = fc.get_kondo_screening(mode="ED")
    assert r["schmidt"] == pytest.approx([0.5]*4, abs=1e-10)
    assert r["M_imp"][0] == pytest.approx(r["M_imp"][1], abs=1e-10)
    assert abs(r["M_env"][0][1]) < 1e-10 # no screening across channels


def test_explicit_impurity_and_bath_match_the_stored_layout():
    fc = wilsonchain.anderson_chain(3, ed=-0.5, U=1.0, tD=0.3, Lambda=4.0)
    a = fc.get_kondo_screening(mode="ED")
    b = fc.get_kondo_screening(imp=0, bath=[1, 2, 3], mode="ED")
    assert a["M_imp"] == pytest.approx(b["M_imp"])
    assert np.allclose(a["M_n"], b["M_n"])


@pytest.mark.parametrize("ed,tD", [(-0.5, 0.3), ([-0.5, -0.5], [0.3, 0.3])])
def test_an_even_chain_warns_since_its_ground_state_is_degenerate(ed, tD):
    with pytest.warns(UserWarning, match="even"):
        wilsonchain.anderson_chain(2, ed=ed, U=1.0, tD=tD, Lambda=4.0)
