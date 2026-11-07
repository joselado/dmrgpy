"""A Kondo impurity of spin S > 1/2 with uniaxial anisotropy on the BCS
Wilson chain (ysr.py): 2S singly occupied sites bound by a Hund coupling,
so that the spin-S multiplet is an exactly decoupled block."""
import fractions

import numpy as np
import pytest

from dmrgpy import ysr


def _spin_observables(fc):
    Sx, Sy, Sz = ysr.impurity_spin(fc)
    return {"S2": Sx*Sx + Sy*Sy + Sz*Sz, "Sz2": Sz*Sz}


def test_spin_parsing_and_layout():
    assert [ysr.two_spin(s) for s in ("1/2", "1", "3/2", "2", 1.5, 2)] == \
        [1, 2, 3, 4, 3, 4]
    for bad in ("1/3", 0, "-1/2"):
        with pytest.raises(ValueError):
            ysr.two_spin(bad)
    fc = ysr.kondo_chain(3, J=0.3, Delta=0.1, S="3/2", Dz=0.02,
                         itensor_version="python")
    assert fc.impurity == [0, 1, 2] and fc.channels == [[3, 4, 5]]
    assert fc.spin == fractions.Fraction(3, 2)


def test_spin_one_half_keeps_its_layout():
    """S=1/2 is the old chain: impurity on site 0, bath from site 1, and
    the tunneling vertices default to that layout."""
    fc = ysr.kondo_chain(3, J=0.6, Delta=0.1, itensor_version="python")
    ref = ysr.kondo_chain(3, J=0.6, Delta=0.1, S=0.5, itensor_version="python")
    assert fc.hamiltonian.op == ref.hamiltonian.op
    assert fc.impurity == [0] and fc.channels == [[1, 2, 3]]
    a = ysr.tunneling_operators(fc, U=0.3)
    b = ysr.tunneling_operators(fc, imp=0, site=1, U=0.3)
    assert all(x.op == y.op for x, y in zip(a, b))
    r = ysr.subgap_energies(fc, mode="ED")
    assert r["E_ysr"] == r["E_odd"] - r["E_even"]
    assert r["E_screened"] == r["E_even"] and set(r["sectors"]) == {0, 1}


@pytest.mark.parametrize("itensor_version", ["python", 3])
@pytest.mark.parametrize("S,Dz,n", [("1", 0.05, 4), ("3/2", -0.04, 4),
                                    ("2", 0.03, 2)])
def test_sector_energies_dmrg_against_ed(S, Dz, n, itensor_version):
    """Every Sz sector, and the impurity in its spin-S multiplet,
    <S_imp^2> = S(S+1), in the solved state of every sector."""
    fc = ysr.kondo_chain(n, J=0.6, Delta=0.1, S=S, Dz=Dz,
                         itensor_version=itensor_version)
    fc.maxm, fc.nsweeps = 64, 15
    obs = _spin_observables(fc)
    ed = ysr.subgap_energies(fc, mode="ED", observables=obs)
    dm = ysr.subgap_energies(fc, observables=obs)
    s = float(fc.spin)
    assert set(ed["sectors"]) == set(range(0, ysr.two_spin(S) + 2))
    for k in ed["sectors"]:
        assert dm["sectors"][k] == pytest.approx(ed["sectors"][k], abs=1e-8)
        for r in (ed, dm):
            assert r["vevs"]["S2"][k].real == pytest.approx(s*(s+1), abs=1e-8)
        assert dm["vevs"]["Sz2"][k].real == \
            pytest.approx(ed["vevs"]["Sz2"][k].real, abs=1e-6)
    assert dm["Sz_ground"] == ed["Sz_ground"]
    assert dm["ground_state"] == ed["ground_state"]


def test_free_spin_limit():
    """At J=0 the spin is free: the sectors of the bath-even parity are
    its levels Dz Sz^2 over the BCS vacuum, exactly."""
    Dz = 0.07
    fc = ysr.kondo_chain(3, J=0.0, Delta=0.1, S="1", Dz=Dz,
                         itensor_version="python")
    r = ysr.subgap_energies(fc, mode="ED")
    assert r["sectors"][2] - r["sectors"][0] == pytest.approx(Dz, abs=1e-12)
    assert r["Sz_ground"] == 0 and r["ground_state"] == "free"
    fc = ysr.kondo_chain(3, J=0.0, Delta=0.1, S="3/2", Dz=Dz,
                         itensor_version="python")
    r = ysr.subgap_energies(fc, mode="ED")
    assert r["sectors"][3] - r["sectors"][1] == pytest.approx(2*Dz, abs=1e-12)
    fc = ysr.kondo_chain(2, J=0.0, Delta=0.1, S="2", Dz=-Dz,
                         itensor_version="python")
    r = ysr.subgap_energies(fc, mode="ED")
    assert set(r["sectors"]) == set(range(6))
    assert r["sectors"][4] - r["sectors"][2] == pytest.approx(-3*Dz, abs=1e-12)
    assert r["sectors"][2] - r["sectors"][0] == pytest.approx(-Dz, abs=1e-12)
    assert r["Sz_ground"] == 4 and r["ground_state"] == "free"


def test_field_picks_the_sign_of_sz():
    """With B the sectors of both signs are scanned, and B S^z with B > 0
    puts the ground state at negative Sz."""
    fc = ysr.kondo_chain(3, J=0.0, Delta=0.1, S="1", B=0.05,
                         itensor_version="python")
    r = ysr.subgap_energies(fc, mode="ED")
    assert set(r["sectors"]) == {-3, -2, -1, 0, 1, 2, 3}
    assert r["Sz_ground"] == -2
    assert r["sectors"][2] - r["sectors"][-2] == pytest.approx(0.1, abs=1e-12)


def test_stm_spectrum_has_the_right_weight_for_spin_one():
    """The integrated spectrum is sum_sigma <{X_sigma, X_sigma^dag}> with
    the total impurity spin in X, less the Lorentzian tails outside the
    window, 2 delta/(pi L) of the weight for poles well inside [-L, L]."""
    fc = ysr.kondo_chain(2, J=0.6, Delta=0.1, S="1", Dz=0.05,
                         itensor_version="python")
    X = ysr.tunneling_operators(fc)
    rule = fc.vev(sum(x*x.get_dagger() + x.get_dagger()*x for x in X),
                  mode="ED").real
    L, delta = 12.0, 0.06
    es = np.linspace(-L, L, 801)
    _, A = ysr.get_stm_spectrum(fc, es, delta=delta, mode="ED",
                                submode="CVM")
    tail = 2*delta/(np.pi*L)
    assert np.trapezoid(A, es) == pytest.approx(rule*(1 - tail), abs=2e-4)


@pytest.mark.parametrize("itensor_version", ["python", 3])
@pytest.mark.parametrize("S,Dz,n", [("1", 0.05, 4), ("1", -0.05, 4),
                                    ("3/2", -0.04, 3)])
def test_stm_spectrum_dmrg_against_ed(S, Dz, n, itensor_version):
    es = np.linspace(-1.2, 1.2, 25)
    fc = ysr.kondo_chain(n, J=0.6, Delta=0.1, S=S, Dz=Dz,
                         itensor_version=itensor_version)
    fc.maxm, fc.nsweeps = 64, 15
    _, ed = ysr.get_stm_spectrum(fc, es, delta=0.08, U=0.3, mode="ED",
                                 submode="CVM")
    _, dm = ysr.get_stm_spectrum(fc, es, delta=0.08, U=0.3)
    assert np.max(np.abs(dm - ed)) < 1e-7
