"""The Kondo impurity on a superconducting substrate (ysr.py): the BCS
Wilson chain, the subgap state, and the STM spectrum against ED."""
import warnings

import numpy as np
import pytest

from dmrgpy import wilsonchain, ysr


def _bdg(eps, t, Delta, v0=0.0):
    """Nambu (f_up, f_dn^dag) matrix of a chain with uniform onsite
    pairing, and a classical exchange field v0 on site 0."""
    n = len(eps)
    h = np.diag(eps) + np.diag(t, 1) + np.diag(t, -1)
    H = np.block([[h, Delta*np.eye(n)], [Delta*np.eye(n), -h]])
    H[0, 0] += v0
    H[n, n] += v0
    return H


def test_pairing_wilson_chain_is_the_bcs_star():
    """Uniform onsite pairing on the Wilson chain gives exactly the
    Nambu Green function of the discretized BCS band at the impurity."""
    xi, w = wilsonchain.star_levels(Lambda=2.0, nint=12)
    eps, t = wilsonchain.wilson_chain(len(xi), Lambda=2.0, nint=12)
    Delta = 0.05
    H = _bdg(eps, t, Delta)
    for om in (0.01 + 0.003j, 0.2 + 0.01j, 0.7 + 0.02j):
        G = np.linalg.inv(om*np.eye(len(H)) - H)[0, 0]
        Gs = np.sum(w*(om + xi)/(om**2 - xi**2 - Delta**2))
        assert abs(G - Gs) < 1e-12*abs(Gs)


@pytest.mark.parametrize("Lambda", [1.5, 2.0, 3.0])
@pytest.mark.parametrize("JS", [0.3, 0.6, 1.0])
def test_classical_ysr_energy(Lambda, JS):
    """A static exchange field J S s^z on the first site binds one subgap
    state at Delta (1-a^2)/(1+a^2), a = pi rho J S/2, where a field v0 on
    the chain acts as v0/A_Lambda in the continuum."""
    Delta = 1e-3
    eps, t = wilsonchain.wilson_chain(80, Lambda=Lambda)
    ev = np.linalg.eigvalsh(_bdg(eps, t, Delta, v0=JS/2))
    a = np.pi*0.5*(JS/2)/ysr.lambda_factor(Lambda)
    pred = Delta*abs(1 - a*a)/(1 + a*a)
    assert abs(np.min(np.abs(ev)) - pred) < 1e-3*Delta


def test_quasiparticle_gap_reaches_delta_on_a_long_chain():
    short = ysr.kondo_chain(4, J=0.2, Delta=1e-2)
    eps, t = wilsonchain.wilson_chain(30)
    assert short.qp_gap > 5e-2
    assert ysr.bath_quasiparticle_gap(eps, t, 1e-2) == pytest.approx(1e-2, rel=1e-4)


def test_the_sign_of_the_gap_does_not_change_reaches_gap():
    """The sign of Delta is a gauge, so the metastable-sector check must
    see a long chain at Delta < 0 as well."""
    for Delta in (1e-2, -1e-2):
        fc = ysr.kondo_chain(12, J=0.2, Delta=Delta, Lambda=4.0)
        assert fc.reaches_gap
        assert fc.qp_gap == pytest.approx(1e-2, rel=1e-2)
    assert not ysr.kondo_chain(10, J=0.2, Delta=-1e-2, Lambda=4.0).reaches_gap


def test_twist_zero_is_the_old_mesh():
    """z=0 keeps the mesh unchanged and every z keeps the weight."""
    xi0, w0 = wilsonchain.star_levels(Lambda=2.0, nint=10)
    pts = 2.0**(-np.arange(11.0))
    assert np.allclose(xi0[:10], (pts[:-1] + pts[1:])/2)
    for z in (0.25, 0.5, 0.75):
        xi, w = wilsonchain.star_levels(Lambda=2.0, nint=10, z=z)
        assert abs(np.sum(w) - 1) < 1e-14 and np.max(xi) < 1.0
    with pytest.raises(ValueError):
        wilsonchain.star_levels(z=1.0)


@pytest.mark.parametrize("itensor_version", ["python", 3])
def test_subgap_energies_dmrg_against_ed(itensor_version):
    fc = ysr.kondo_chain(4, J=0.6, Delta=0.1, itensor_version=itensor_version)
    fc.maxm, fc.nsweeps = 60, 15
    ed = ysr.subgap_energies(fc, mode="ED")
    dm = ysr.subgap_energies(fc)
    for k in ("E_even", "E_odd"):
        assert dm[k] == pytest.approx(ed[k], abs=1e-8)


def test_gap_drives_the_singlet_doublet_transition():
    """Small gap: the spin is screened. Large gap: it is free."""
    small = ysr.subgap_energies(ysr.kondo_chain(5, J=0.6, Delta=0.01),
                                mode="ED")
    large = ysr.subgap_energies(ysr.kondo_chain(5, J=0.6, Delta=1.0),
                                mode="ED")
    assert small["ground_state"] == "singlet"
    assert large["ground_state"] == "doublet"


def test_stm_spectrum_has_the_right_weight():
    """The integrated spectrum is sum_sigma <{X_sigma, X_sigma^dag}>,
    which at U=0 is 3/8 - <S.s_1>: the tau matrices do not commute, so it
    is not a constant."""
    fc = ysr.kondo_chain(3, J=0.6, Delta=0.1, itensor_version="python")
    X = ysr.tunneling_operators(fc)
    rule = fc.vev(sum(x*x.get_dagger() + x.get_dagger()*x for x in X),
                  mode="ED").real
    ss = fc.Sx[0]*fc.Sx[1] + fc.Sy[0]*fc.Sy[1] + fc.Sz[0]*fc.Sz[1]
    assert rule == pytest.approx(3/8 - fc.vev(ss, mode="ED").real, abs=1e-12)
    es = np.linspace(-30, 30, 6001)
    _, A = ysr.get_stm_spectrum(fc, es, delta=0.02, mode="ED",
                                submode="CVM")
    assert np.trapezoid(A, es) == pytest.approx(rule, abs=1e-3)


@pytest.mark.parametrize("itensor_version", ["python", 3])
@pytest.mark.parametrize("U", [0.0, 0.3])
def test_stm_spectrum_dmrg_against_ed(itensor_version, U):
    es = np.linspace(-1.2, 1.2, 25)
    fc = ysr.kondo_chain(4, J=0.6, Delta=0.1, itensor_version=itensor_version)
    fc.maxm, fc.nsweeps = 60, 15
    _, ed = ysr.get_stm_spectrum(fc, es, delta=0.08, U=U, mode="ED",
                                 submode="CVM")
    _, dm = ysr.get_stm_spectrum(fc, es, delta=0.08, U=U)
    assert np.max(np.abs(dm - ed)) < 1e-8


def test_tip_convolution_normal_tip_is_identity():
    es = np.linspace(-1, 1, 801)
    A = 1 + 0.3*np.exp(-(es - 0.2)**2/0.01)
    Vs, didv = ysr.tip_convolution(es, A, 0.0, Vs=es[100:-100])
    assert np.max(np.abs(didv - A[100:-100])) < 1e-10


@pytest.mark.parametrize("Delta,ratio", [(0.003, 0.5606), (0.03, -0.2302)])
def test_kondo_ysr_transition_on_a_long_chain(Delta, ratio):
    """rho J = 0.2 in the continuum on 14 sites at Lambda=4: screened at
    Delta=0.003, free at Delta=0.03. On d22, 20 sites at maxm=128 gave
    Delta_c/T_K = 4.4 to 5.0 for rho J = 0.09 to 0.19, extrapolating to
    3.74 at rho J -> 0 (NRG: 3.7)."""
    Lambda = 4.0
    fc = ysr.kondo_chain(14, J=0.4*ysr.lambda_factor(Lambda), Delta=Delta,
                         Lambda=Lambda, itensor_version=3)
    fc.maxm, fc.nsweeps = 64, 20
    assert fc.reaches_gap
    r = ysr.subgap_energies(fc)
    assert r["E_ysr"]/Delta == pytest.approx(ratio, abs=1e-3)
