"""The ground state of definite fermion parity for the STM spectrum
(ysrtk/paritystate.py): solved inside its 2Sz sector, promoted to dense
indices and handed to the TASK correlator. Every tunneling vertex changes
the parity, so on this state the elastic weight vanishes exactly, while
the unrestricted solve near the singlet-doublet transition can stop at a
mixture of the two parities."""
import warnings

import numpy as np
import pytest

from dmrgpy import ysr
from dmrgpy.ysrtk import paritystate, surrogate, ysrcurrent, zaverage


@pytest.mark.parametrize("itensor_version", ["python", 3])
@pytest.mark.parametrize("U", [0.0, 0.3])
def test_task_on_the_promoted_parity_state_against_ed(itensor_version, U):
    """The chain of test_ysr.py::test_stm_spectrum_dmrg_against_ed. Run on
    d22: 2.6e-11 (v3) and 8.2e-10 ("python") against ED."""
    es = np.linspace(-1.2, 1.2, 25)
    fc = ysr.kondo_chain(4, J=0.6, Delta=0.1, itensor_version=itensor_version)
    fc.maxm, fc.nsweeps = 60, 15
    _, ed = ysr.get_stm_spectrum(fc, es, delta=0.08, U=U, mode="ED",
                                 submode="CVM")
    c, info = paritystate.parity_ground_state(fc)
    assert c is not fc and c.conserved_sector is None
    assert fc.conserved_sector is None # fc itself is untouched
    assert info["Sz_ground"] == 1 and info["parity"] == "odd"
    assert info["E0"] == pytest.approx(min(info["sectors"].values()))
    _, dm = ysr.get_stm_spectrum(c, es, delta=0.08, U=U)
    assert np.max(np.abs(dm - ed)) < 1e-8
    p = zaverage.stm_poles(c, U=U, nkry=100)
    assert p["elastic"] < 1e-12 and p["missing_weight"] < 1e-12


_NEAR = dict(J=0.6*ysr.lambda_factor(4.0), Delta=0.1, Lambda=4.0)
_EPS_ED = 0.019031164828 # E_ysr/Delta by ED, sector solves on the same chain


def test_parity_state_removes_the_mixture_near_the_transition():
    """8 sites at Lambda=4 near the singlet-doublet transition, the chain of
    the 'Width of a YSR peak' note. At maxm=60, nsweeps=15 the plain v3
    solve stopped at a mixture in 4 runs of 4 on d22 (elastic weight 1.2e-2,
    9.7e-5, 3.8e-3, 3.6e-6 of the total). In the parity sector the elastic
    weight is exactly 0 and rate_parameters gives the converged u2, v2."""
    Delta = _NEAR["Delta"]
    fc = ysr.kondo_chain(8, itensor_version=3, **_NEAR)
    fc.maxm, fc.nsweeps = 60, 15
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        p = zaverage.stm_poles(fc, U=0.3, nkry=100, parity=True)
        r = ysrcurrent.rate_parameters(p, Delta, Lambda=4.0)
    assert [str(x.message) for x in w
            if issubclass(x.category, RuntimeWarning)] == []
    assert p["elastic"] < 1e-10*np.sum(np.abs(p["weight"]))
    assert p["parity"]["Sz_ground"] == 0
    assert p["parity"]["E_ysr"]/Delta == pytest.approx(_EPS_ED, abs=1e-6)
    om, wt = p["omega"], p["weight"]
    sub = (np.abs(om) < Delta) & (np.abs(wt) > 1e-6*np.max(np.abs(wt)))
    for sign in (1.0, -1.0):
        got = om[sub & (sign*om > 0.0)]
        assert np.allclose(sign*got/Delta, _EPS_ED, atol=1e-6)
    assert r["eps"]/Delta == pytest.approx(_EPS_ED, abs=1e-6)
    assert r["u2"] == pytest.approx(0.333787, abs=1e-5)
    assert r["v2"] == pytest.approx(0.196261, abs=1e-5)


@pytest.mark.parametrize("itensor_version", ["python", 3])
@pytest.mark.parametrize("S,Dz,n,cutoff,tol", [("1", 0.05, 4, None, 1e-7),
                                               ("1", -0.05, 4, 0.0, 1e-6),
                                               ("3/2", -0.04, 3, 0.0, 1e-8)])
def test_higher_spin_parity_state_against_ed(S, Dz, n, cutoff, tol,
                                             itensor_version):
    """Spin S as 2S sites with anisotropy: the sectors up to |2Sz| = 2S+1,
    the seeded v3 solve and the charge check, then TASK on the promoted
    state against ED. With Dz < 0 the ground state is at 2Sz = +-2, and
    the sector state at the chain's cutoff 1e-12 has a tangent space of
    562 (S=1) and 205 (S=3/2) against 786 and 297 for the plain solve, which
    mixes the two members, so its TASK spectrum is off ED by 8.7e-4 and
    4.8e-2 of the peak (d22, both backends; the subgap pole is exact to
    4e-10). With cutoff=0 the tangent space is the whole Hilbert space
    (4095 on "python", 2870 on v3) and the spectrum is ED's: 2.6e-7 and
    1.6e-7 for S=1, 8e-11 and 2e-10 for S=3/2."""
    es = np.linspace(-1.2, 1.2, 25)
    fc = ysr.kondo_chain(n, J=0.6, Delta=0.1, S=S, Dz=Dz,
                         itensor_version=itensor_version)
    fc.maxm, fc.nsweeps = 64, 15
    if cutoff is not None:
        fc.cutoff = cutoff
    _, ed = ysr.get_stm_spectrum(fc, es, delta=0.08, U=0.3, mode="ED",
                                 submode="CVM")
    ref = ysr.subgap_energies(fc, mode="ED")
    c, info = paritystate.parity_ground_state(fc)
    assert sorted(info["sectors"]) == sorted(ref["sectors"])
    assert info["Sz_ground"] == ref["Sz_ground"]
    assert info["E_ysr"] == pytest.approx(ref["E_ysr"], abs=1e-8)
    _, dm = ysr.get_stm_spectrum(c, es, delta=0.08, U=0.3)
    assert np.max(np.abs(dm - ed)) < tol


def test_ed_mode_returns_a_plain_copy_with_the_sector_energies():
    """On mode="ED" the full diagonalization already has a definite parity,
    and an ED chain with a sector set would restrict the vertex to it."""
    fc = ysr.kondo_chain(3, J=0.6, Delta=0.1, itensor_version="python")
    c, info = paritystate.parity_ground_state(fc, mode="ED")
    assert c is not fc and c.conserved_sector is None
    ref = ysr.subgap_energies(fc, mode="ED")
    assert info["E_ysr"] == pytest.approx(ref["E_ysr"], abs=1e-12)
    assert info["Sz_ground"] == ref["Sz_ground"]


def test_parity_state_needs_a_gap_and_both_parities():
    fc = ysr.kondo_chain(3, J=0.6, Delta=0.0, itensor_version="python")
    with pytest.raises(ValueError, match="gap"):
        paritystate.parity_ground_state(fc)
    fc = ysr.kondo_chain(3, J=0.6, Delta=0.1, itensor_version="python")
    with pytest.raises(ValueError, match="both parities"):
        paritystate.parity_ground_state(fc, sectors=[0, 2])


def test_default_sectors_follow_subgap_energies_and_sector_energies():
    py = dict(itensor_version="python")
    assert paritystate.default_sectors(ysr.kondo_chain(2, Delta=0.1, **py)) \
        == [0, 1]
    assert paritystate.default_sectors(ysr.kondo_chain(
        2, Delta=0.1, S="1", Dz=0.05, **py)) == [0, 1, 2, 3]
    assert paritystate.default_sectors(ysr.kondo_chain(
        2, Delta=0.1, B=0.01, **py)) == [-2, -1, 0, 1, 2]
    assert paritystate.default_sectors(surrogate.surrogate_kondo_chain(
        2, J=0.6, Delta=0.1, L=1, B=0.01, **py)) == list(range(-3, 4))


def test_stm_poles_warns_when_the_tangent_space_misses_the_vertex():
    """J = 0: a product ground state, and the composite vertex leaves the
    one-site tangent space, so TASK finds none of the weight."""
    fc = ysr.kondo_chain(3, J=0.0, Delta=0.1, itensor_version="python")
    fc.maxm, fc.nsweeps = 32, 10
    with pytest.warns(RuntimeWarning, match="bdg_chain_dos"):
        p = zaverage.stm_poles(fc, nkry=50)
    assert p["missing_weight"] > 0.5*(p["missing_weight"]
                                      + np.sum(np.abs(p["weight"])))


# A surrogate dimer next to the crossing of the screened and the free spin:
# the lowest state of odd parity is 1.709e-3 above the ground state (ED).
_DIMER = dict(J=0.671, Delta=0.1, L=3)
_DIMER_ES = np.array([-0.05, -0.0017088, 0.0, 0.0017088, 0.05])


@pytest.fixture(scope="module")
def dimer_ed():
    fc = surrogate.surrogate_kondo_chain(2, **_DIMER)
    _, A = surrogate.get_stm_spectrum(fc, _DIMER_ES, spin=0, delta=0.005,
                                      U=0.3, mode="ED", submode="CVM")
    return A


@pytest.mark.parametrize("itensor_version", ["python", 3])
def test_surrogate_dimer_near_a_parity_degeneracy(itensor_version, dimer_ed):
    """J12 = 0, rho J = 0.34: one spin screened against both, odd gap
    1.709e-3 by ED. On d22 three plain v3 solves were 9.4e-5 to 5.2e-4 above
    the ED energy and off by 1.1e-2 to 4.0e-2 of the peak (16.15); the
    parity state was off by 3.5e-12 of it."""
    ed = dimer_ed
    fc = surrogate.surrogate_kondo_chain(2, itensor_version=itensor_version,
                                         **_DIMER)
    fc.maxm, fc.nsweeps = 64, 20
    _, A = surrogate.get_stm_spectrum(fc, _DIMER_ES, spin=0, delta=0.005,
                                      U=0.3, parity=True)
    assert np.max(np.abs(A - ed)) < 1e-6*np.max(ed)
    assert fc.conserved_sector is None


# The same dimer with the two spins coupled, by an exchange and by a hopping
# between the first sites of their baths. Its ED (CVM) spectrum at _DIMER_ES,
# delta=0.005, U=0.3, computed on d22 (80 s, too slow for a fixture).
_COUPLED = dict(_DIMER, J12=0.01, t_b=0.01)
_COUPLED_ED = np.array([0.13083578, 6.37471529, 5.97715881, 6.19345192,
                        0.13088216])


@pytest.mark.parametrize("itensor_version,nsweeps", [("python", 40), (3, 20)])
def test_coupled_dimer_parity_state(itensor_version, nsweeps):
    """J12 = t_b = 0.01. On "python" the sector solve converges more slowly
    than the plain one here, so it gets twice the sweeps (next test)."""
    ed = _COUPLED_ED
    fc = surrogate.surrogate_kondo_chain(2, itensor_version=itensor_version,
                                         **_COUPLED)
    fc.maxm, fc.nsweeps = 64, nsweeps
    _, A = surrogate.get_stm_spectrum(fc, _DIMER_ES, spin=0, delta=0.005,
                                      U=0.3, parity=True)
    assert np.max(np.abs(A - ed)) < 1e-5*np.max(ed)


@pytest.mark.xfail(strict=True, reason="the \"python\" sector solve of the "
                   "coupled dimer is not converged at nsweeps=20")
def test_coupled_dimer_parity_state_python_at_twenty_sweeps():
    """On d22 the parity route was off ED by 1.1e-4 of the peak here, with
    its energy within 1e-10 of the ED sector energy, while the plain route
    at the same settings was 5e-7 off. The state stays in its sector
    (<Sz^2> 6e-15) but its sum rule <X^dag X> is 1.3e-5 off, so the error
    is an admixture of another state of the same sector, not a leak out
    of it."""
    ed = _COUPLED_ED
    fc = surrogate.surrogate_kondo_chain(2, itensor_version="python",
                                         **_COUPLED)
    fc.maxm, fc.nsweeps = 64, 20
    _, A = surrogate.get_stm_spectrum(fc, _DIMER_ES, spin=0, delta=0.005,
                                      U=0.3, parity=True)
    assert np.max(np.abs(A - ed)) < 1e-5*np.max(ed)


def test_zaverage_spectrum_with_parity_keeps_the_subgap_pole():
    """The 5-site chain of test_ysr_zaverage.py's short-chain test: with
    parity=True the subgap pole is the same, the elastic weight is 0 and
    info carries the sector energies of each mesh."""
    es = np.linspace(-0.2, 0.2, 101)
    with pytest.warns(RuntimeWarning, match="does not reach"):
        _, A, info = zaverage.zaverage_stm_spectrum(5, es, J=0.15,
                Delta=0.05, Lambda=4.0, nz=1, b=0.4, U=0.3, delta_sub=1e-3,
                maxm=64, nsweeps=20, nkry=100, parity=True, parity_gap=True,
                itensor_version="python")
    assert sorted(info["parity"][0]) == [0, 1]
    assert info["elastic"][0] < 1e-12
    sub = info["subgap"][0]
    big = np.abs(sub["weight"]) > 1e-8
    assert np.allclose(np.abs(sub["omega"][big]), 0.04189, atol=1e-4)
    assert np.allclose(np.abs(sub["omega"][big]), abs(info["E_ysr"][0]),
                       atol=1e-8)
    assert info["E0"][0] == pytest.approx(min(info["parity"][0].values()))
