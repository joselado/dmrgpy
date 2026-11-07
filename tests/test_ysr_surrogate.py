"""Surrogate BCS baths (ysrtk/surrogate.py): the fit against the published
parameters, the chain map, the Kondo map on the classical YSR energy, and
one and two spins against ED."""
import numpy as np
import pytest
from scipy.optimize import brentq

from dmrgpy import ysr
from dmrgpy.ysrtk import surrogate as su

# gamma and xi (units of Delta) of github.com/virgilb91/surrogate_models,
# files gamma_L=*_B=10_wc=10.dat and xi_L=*_B=10_wc=10.dat (MIT license),
# the fit of arXiv:2307.11646 Table I at D = wc = 10 Delta
PUBLISHED = {
    1: ([0.9617925826975624], [0.0]),
    2: ([1.2641891616930312]*2, [-1.309892193517106, 1.309892193517106]),
    3: ([0.5079990272016913, 1.8454385496053414, 1.8454385496053414],
        [0.0, 2.7545651912448323, -2.7545651912448323]),
    4: ([0.5258044702613829]*2 + [2.113958393542536]*2,
        [0.7066985291471469, -0.7066985291471469, 4.122688551050398,
         -4.122688551050398]),
}


@pytest.mark.parametrize("L", [1, 2, 3, 4])
def test_fit_gives_the_published_parameters(L):
    f = su.fit_surrogate(L, Delta=1.0, D=10.0)
    g0, x0 = map(np.array, PUBLISHED[L])
    o, o0 = np.argsort(f["xi"]), np.argsort(x0)
    assert np.allclose(f["xi"][o], x0[o0], atol=1e-6)
    assert np.allclose(f["gamma"][o], g0[o0], atol=1e-6)


def test_fit_scales_with_the_gap():
    a = su.fit_surrogate(3, Delta=1.0, D=10.0)
    b = su.fit_surrogate(3, Delta=0.05, D=0.5)
    assert np.allclose(b["xi"], 0.05*a["xi"]) and np.allclose(b["gamma"], 0.05*a["gamma"])
    assert b["W"] == pytest.approx(a["W"])


def test_chain_is_the_surrogate_star():
    """The Lanczos chain with uniform pairing has the Nambu Green function
    of the star at its first site."""
    Delta = 0.1
    f = su.fit_surrogate(4, Delta=Delta, D=1.0)
    eps, t = su.star_to_chain(f["xi"], f["gamma"])
    assert np.all(eps == 0.0)
    n = len(eps)
    h = np.diag(eps) + np.diag(t, 1) + np.diag(t, -1)
    H = np.block([[h, Delta*np.eye(n)], [Delta*np.eye(n), -h]])
    w = f["gamma"]/np.sum(f["gamma"])
    for om in (0.01 + 0.003j, 0.2 + 0.01j, 0.7 + 0.02j):
        G = np.linalg.inv(om*np.eye(2*n) - H)[0, 0]
        Gs = np.sum(w*(om + f["xi"])/(om**2 - f["xi"]**2 - Delta**2))
        assert abs(G - Gs) < 1e-12*abs(Gs)


def _classical_ysr(v0, Delta, D):
    """Bound state of a static field v0 on the local orbital of a flat BCS
    band of half-width D: 1 = v0 g(E) (E + Delta)."""
    rho = 1/(2*D)
    def f(E):
        s = np.sqrt(Delta**2 - E**2)
        return 1 - v0*2*rho*np.arctan(D/s)/s*(E + Delta)
    return abs(brentq(f, -Delta*(1 - 1e-13), Delta*(1 - 1e-13)))


@pytest.mark.parametrize("DD", [10.0, 1000.0])
def test_classical_ysr_energy_fixes_the_weight(DD):
    """A field v on the local orbital acts as v W on the normalized first
    site. With that, the classical YSR energy is right at L=4 to a few
    percent of Delta, at any D/Delta."""
    D, v0 = 1.0, 1.0
    Delta = D/DD
    f = su.fit_surrogate(4, Delta=Delta, D=D)
    eps, t = su.star_to_chain(f["xi"], f["gamma"])
    n = len(eps)
    h = np.diag(eps) + np.diag(t, 1) + np.diag(t, -1)
    H = np.block([[h, Delta*np.eye(n)], [Delta*np.eye(n), -h]])
    H[0, 0] += v0*f["W"]
    H[n, n] += v0*f["W"]
    e = np.min(np.abs(np.linalg.eigvalsh(H)))
    assert e == pytest.approx(_classical_ysr(v0, Delta, D), abs=0.02*Delta)


def test_one_spin_matches_the_ysr_convention():
    fc = su.surrogate_kondo_chain(1, J=0.6, Delta=0.1, L=3)
    r = su.sector_energies(fc, mode="ED")
    assert r["bath_parity"] == {0: "odd", 1: "even"}
    e = ysr.subgap_energies(fc, mode="ED")
    assert e["E_ysr"] == pytest.approx(r["sectors"][1][0] - r["sectors"][0][0],
                                       abs=1e-12)


def test_field_solves_both_signs_of_sz():
    """With B != 0 the lowest member of the doublet is at 2Sz = -1, which
    the default sectors must contain (0.18 Delta off with (0, 1) alone)."""
    fc = su.surrogate_kondo_chain(1, J=0.6, Delta=0.1, L=3, B=0.02)
    r = su.sector_energies(fc, mode="ED")
    e = ysr.subgap_energies(fc, mode="ED")
    assert -1 in r["sectors"] and r["Sz_ground"] == e["Sz_ground"] == -1
    assert r["E_odd"] == pytest.approx(e["E_odd"], abs=1e-12)
    # pinned on one machine; another BLAS moves the converged E0 by 1.4e-9
    assert r["E0"] == pytest.approx(-0.7289046253, abs=1e-8)


def test_degenerate_multiplet_takes_the_smallest_sz():
    fc = su.surrogate_kondo_chain(2, J=0.3, Delta=0.1, L=2, J12=-0.5)
    r = su.sector_energies(fc, mode="ED", sectors=(0, 1, 2))
    assert r["sectors"][2][0] == pytest.approx(r["sectors"][0][0], abs=1e-12)
    assert r["Sz_ground"] == 0


def test_a_fit_made_elsewhere_is_refused():
    f = su.fit_surrogate(3, Delta=0.1, D=1.0)
    for kw in (dict(Delta=0.2), dict(Delta=0.1, L=4), dict(Delta=0.1, D=2.0)):
        with pytest.raises(ValueError):
            su.surrogate_kondo_chain(1, fit=f, **kw)
    fc = su.surrogate_kondo_chain(1, fit=f, Delta=0.1)
    assert fc.ns == 4


def test_decoupled_dimer_is_two_monomers():
    mono = su.sector_energies(su.surrogate_kondo_chain(1, J=0.6, Delta=0.1, L=2),
                              mode="ED")["sectors"]
    di = su.sector_energies(su.surrogate_kondo_chain(2, J=0.6, Delta=0.1, L=2),
                            mode="ED", sectors=(0, 1, 2))["sectors"]
    Es, Ed = mono[0][0], mono[1][0]
    assert di[0][0] == pytest.approx(min(2*Es, 2*Ed), abs=1e-10)
    assert di[1][0] == pytest.approx(Es + Ed, abs=1e-10)
    assert di[2][0] == pytest.approx(2*Ed, abs=1e-10)


@pytest.mark.parametrize("itensor_version", ["python", 3])
def test_dimer_dmrg_against_ed(itensor_version):
    fc = su.surrogate_kondo_chain(2, J=1.5, V=1.0, Delta=0.1, L=2, t_b=0.05,
                                  J12=0.03, itensor_version=itensor_version)
    fc.maxm, fc.nsweeps = 60, 15
    ed = su.sector_energies(fc, mode="ED", sectors=(0, 1, 2))
    dm = su.sector_energies(fc, sectors=(0, 1, 2))
    for s in (0, 1, 2):
        assert dm["sectors"][s][0] == pytest.approx(ed["sectors"][s][0], abs=1e-8)


def test_dimer_excited_states_v3_against_ed():
    fc = su.surrogate_kondo_chain(2, J=1.5, V=1.0, Delta=0.1, L=2, t_b=0.05,
                                  itensor_version=3)
    fc.maxm, fc.nsweeps = 60, 15
    ed = su.sector_energies(fc, mode="ED", sectors=(1,), nstates=2)
    dm = su.sector_energies(fc, sectors=(1,), nstates=2)
    assert np.allclose(dm["sectors"][1], ed["sectors"][1], atol=1e-7)


def test_python_excited_states_in_a_sector_refused():
    fc = su.surrogate_kondo_chain(2, J=1.5, Delta=0.1, L=2,
                                  itensor_version="python")
    with pytest.raises(NotImplementedError):
        su.sector_energies(fc, sectors=(1,), nstates=2)


def test_ysr_states_hybridize_only_without_particle_hole_symmetry():
    """The two odd states of a dimer split linearly in t_b at V != 0 and
    not at all at particle-hole symmetry, V = 0 (in the zero-bandwidth
    model of arXiv:2203.15011, eq. A8, the splitting is t_b |u^2 - v^2|)."""
    def split(V, tb):
        fc = su.surrogate_kondo_chain(2, J=1.5, V=V, Delta=0.1, L=2, t_b=tb)
        e = su.sector_energies(fc, mode="ED", sectors=(1,), nstates=2)["sectors"][1]
        return e[1] - e[0]
    assert split(0.0, 0.04) < 1e-10
    s1, s2 = split(1.0, 0.01), split(1.0, 0.02)
    assert s1 > 1e-3 and s2/s1 == pytest.approx(2.0, rel=2e-2)


def test_strong_antiferromagnetic_exchange_empties_the_gap():
    """Two spins bound in a singlet have no subgap state: the lowest odd
    excitation goes to the gap edge of the bath (arXiv:1008.1734 VII)."""
    def exc(J12):
        fc = su.surrogate_kondo_chain(2, J=0.9, Delta=0.1, L=3, J12=J12)
        r = su.sector_energies(fc, mode="ED")
        return (r["sectors"][1][0] - r["sectors"][0][0])/fc.qp_gap
    assert exc(0.0) == pytest.approx(0.4428, abs=1e-3)
    assert exc(1.0) == pytest.approx(0.9444, abs=1e-3)


_STM = dict(J=2.0, V=1.0, Delta=0.1, L=2, t_b=0.02)
_ES = np.linspace(-0.2, 0.2, 21)


@pytest.fixture(scope="module")
def stm_ed():
    fc = su.surrogate_kondo_chain(2, **_STM)
    return su.get_stm_spectrum(fc, _ES, spin=0, delta=0.01, U=0.3,
                               mode="ED", submode="CVM")[1]


@pytest.mark.parametrize("itensor_version", ["python", 3])
def test_dimer_stm_spectrum_dmrg_against_ed(itensor_version, stm_ed):
    """Gapped ground state (first excitation 0.65 Delta): TASK is exact."""
    fc = su.surrogate_kondo_chain(2, itensor_version=itensor_version, **_STM)
    fc.maxm, fc.nsweeps = 64, 20
    _, A = su.get_stm_spectrum(fc, _ES, spin=0, delta=0.01, U=0.3)
    assert np.max(np.abs(A - stm_ed)) < 1e-6


_FOUR = dict(J=1.5, V=1.0, Delta=0.1, L=1, t_b=0.05)


@pytest.fixture(scope="module")
def four_ed():
    return su.sector_energies(su.surrogate_kondo_chain(4, **_FOUR), mode="ED")


@pytest.mark.parametrize("itensor_version", ["python", 3])
def test_four_spins_sector_solves_hold_one_electron_per_spin(itensor_version,
                                                             four_ed):
    """H conserves the charge of every spin site, so a sector solve must
    start with one electron on each: from a random v3 sector start both
    sectors stopped at E = +4.40, with one spin site empty or doubly
    occupied. v3 converges slowly in the odd sector, where one
    quasiparticle has to spread over the chain: 1.2e-3 above ED after 20
    sweeps, exact after 80, so there it is only checked as an upper bound
    with the right charges (sector_energies raises otherwise)."""
    fc = su.surrogate_kondo_chain(4, itensor_version=itensor_version, **_FOUR)
    fc.maxm, fc.nsweeps = 64, 20
    dm = su.sector_energies(fc)
    ed = four_ed["sectors"]
    assert dm["sectors"][0][0] == pytest.approx(ed[0][0], abs=1e-8)
    if itensor_version == "python":
        assert dm["sectors"][1][0] == pytest.approx(ed[1][0], abs=1e-8)
    else:
        assert ed[1][0] - 1e-10 < dm["sectors"][1][0] < ed[1][0] + 3e-3
