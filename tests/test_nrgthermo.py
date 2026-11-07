"""Tests of nrgthermo: the NRG thermodynamics of the Kondo impurity on the
Wilson chain of ysr.kondo_chain, and the Kondo temperature in Wilson's
definition."""

import numpy as np
import pytest

from dmrgpy import nrgthermo, ysr, wilsonchain

A4 = ysr.lambda_factor(4.0)


def dense_kondo_spectrum(n, J, Lambda, B=0.0, z=0.0, mesh="mean"):
    """Full spectrum of the spin-1/2 Kondo impurity (a bare doublet) on
    the Wilson chain of n sites, by dense ED with Jordan-Wigner modes."""
    eps, t = wilsonchain.wilson_chain(n, Lambda=Lambda, band=(-1.0, 1.0),
                                      z=z, mesh=mesh)
    nm = 2*n # modes (site, spin), up before dn
    Z, a = np.diag([1.0, -1.0]), np.array([[0.0, 1.0], [0.0, 0.0]])
    c = []
    for j in range(nm):
        m = np.eye(2) # the impurity spin first
        for o in [Z]*j + [a] + [np.eye(2)]*(nm-j-1):
            m = np.kron(m, o)
        c.append(m)
    one = np.eye(2**nm)
    S = [np.kron(o, one) for o in (np.array([[0, 0.5], [0.5, 0]]),
                                    np.array([[0, -0.5j], [0.5j, 0]]),
                                    np.diag([0.5, -0.5]))]
    up, dn = c[0], c[1]
    s0 = [0.5*(up.T@dn + dn.T@up), -0.5j*(up.T@dn - dn.T@up),
          0.5*(up.T@up - dn.T@dn)]
    H = J*sum(S[i]@s0[i] for i in range(3)) + B*S[2]
    for k in range(n):
        H = H + eps[k]*(c[2*k].T@c[2*k] + c[2*k+1].T@c[2*k+1])
    for k in range(n-1):
        for s in (0, 1):
            h = c[2*k+s].T@c[2*k+2+s]
            H = H + t[k]*(h + h.T)
    return np.linalg.eigvalsh(H)


def test_untruncated_ground_energy_equals_kondo_chain_ed():
    J = 0.2*A4
    r = nrgthermo.kondo_nrg(3, J=J, Lambda=4.0, keep=None)
    fc = ysr.kondo_chain(3, J=J, Delta=0.0, Lambda=4.0,
                         itensor_version="python")
    assert r.E0[-1] == pytest.approx(fc.gs_energy(mode="ED"), abs=1e-12)


@pytest.mark.parametrize("B,z,mesh,J", [(0.0, 0.0, "mean", 0.2*A4),
                                        (0.01, 0.3, "mean", 0.2*A4),
                                        (0.0, 0.5, "zp", 0.2)])
def test_untruncated_spectrum_equals_dense_ed(B, z, mesh, J):
    r = nrgthermo.kondo_nrg(3, J=J, Lambda=4.0, keep=None, B=B, z=z,
                            mesh=mesh)
    ed = dense_kondo_spectrum(3, J, 4.0, B=B, z=z, mesh=mesh)
    assert len(r.spectrum) == 128
    assert np.max(np.abs(np.sort(r.spectrum) - ed)) < 1e-12


@pytest.mark.parametrize("z,mesh", [(0.0, "mean"), (0.3, "mean"), (0.5, "zp")])
def test_free_impurity_untruncated(z, mesh):
    """At J = 0 the impurity is a free moment: T chi_imp = 1/4 (Curie)
    and S_imp = ln 2 at every shell, also with onsite energies (z > 0),
    which the free-chain subtraction must remove."""
    r = nrgthermo.kondo_nrg(5, J=0.0, Lambda=4.0, keep=None, z=z, mesh=mesh)
    assert np.max(np.abs(r.Tchi_imp - 0.25)) < 1e-12
    assert np.max(np.abs(r.S_imp - np.log(2.0))) < 1e-12


def test_free_impurity_truncated():
    """The same with truncation, keep=300: the states lost are Boltzmann
    suppressed (measured 1.8e-7 and 1.5e-5)."""
    r = nrgthermo.kondo_nrg(30, J=0.0, Lambda=4.0, keep=300)
    assert np.max(np.abs(r.Tchi_imp - 0.25)) < 1e-6
    assert np.max(np.abs(r.S_imp - np.log(2.0))) < 1e-4


def test_kondo_temperature_wilson_regression():
    """rho J = 0.1 at Lambda = 4, mean mesh, z = 0: the impurity is
    screened (S_imp -> 0) and T_K(Wilson) = 1.33918e-5 (keep 400, 800 and
    1600 agree to 1e-6, beta_bar 0.5 and 1.0 move it by 1.4e-3)."""
    r = nrgthermo.kondo_nrg(40, J=0.2*A4, Lambda=4.0, keep=400)
    k = nrgthermo.kondo_temperature(r)
    assert k["T_K"] == pytest.approx(1.33918e-5, rel=1e-4)
    assert abs(r.S_imp[0] - np.log(2.0)) < 0.01
    assert abs(r.S_imp[-1]) < 0.05 # an even-odd oscillation around 0
    assert abs(r.S_imp[-1] + r.S_imp[-2]) < 1e-3


def test_local_susceptibility_against_dmrg():
    """chi_loc(T -> 0) with B = 1.2e-6 on the impurity, against the DMRG
    value 2(E0 - E(B))/B^2 = 8604 of kondo_chain(27), Nf = 28, Sz = 0,
    maxm 128 (T_K = 1.19949e-5). NRG gives 8605.4."""
    r = nrgthermo.kondo_nrg(40, J=0.2*A4, Lambda=4.0, keep=400, B=1.2e-6)
    k = nrgthermo.kondo_temperature(r, definition="local")
    assert k["chi0"] == pytest.approx(8604.0, rel=1e-3)
    assert k["T_K"] == pytest.approx(1.19949e-5, rel=1e-3)


def test_unknown_definition_raises():
    r = nrgthermo.kondo_nrg(10, J=0.4*A4, Lambda=4.0, keep=100)
    with pytest.raises(ValueError):
        nrgthermo.kondo_temperature(r, definition="Wilson")
    with pytest.raises(ValueError):
        nrgthermo.kondo_temperature(r, definition="local") # B = 0


def test_free_chain_levels_against_dense():
    """The free-chain levels grown site by site give the thermodynamics of
    a dense diagonalization where that one is still accurate."""
    from dmrgpy.nrgthermo import _free_chain_levels, _free_thermo
    for L, n, mesh, z in [(4.0, 25, "mean", 0.0), (2.0, 40, "zp", 0.5)]:
        eps, t = wilsonchain.wilson_chain(n, Lambda=L, band=(-1.0, 1.0),
                                          z=z, mesh=mesh)
        lev = _free_chain_levels(eps, t)
        for N in range(1, n):
            beta = 0.75/abs(t[N-1])
            h = np.diag(eps[:N+1]) + np.diag(t[:N], 1) + np.diag(t[:N], -1)
            a = _free_thermo(lev[N], beta)
            b = _free_thermo(np.linalg.eigvalsh(h), beta)
            assert a == pytest.approx(b, abs=1e-9)


def test_long_chain_keeps_free_moment_and_kondo_temperature():
    """Down to T_N of 1e-21 D: at J = 0 the moment stays free, and T_K
    does not depend on the length of the chain."""
    r0 = nrgthermo.kondo_nrg(70, J=0.0, Lambda=4.0, keep=400)
    assert np.max(np.abs(r0.Tchi_imp - 0.25)) < 1e-6
    assert np.max(np.abs(r0.S_imp - np.log(2.0))) < 1e-4
    r = nrgthermo.kondo_nrg(60, J=0.2*A4, Lambda=4.0, keep=400)
    assert nrgthermo.kondo_temperature(r)["T_K"] == pytest.approx(1.33918e-5,
                                                                  rel=1e-4)


def test_crossing_definition():
    r = nrgthermo.kondo_nrg(40, J=0.2*A4, Lambda=4.0, keep=400)
    kw = nrgthermo.kondo_temperature(r)
    kc = nrgthermo.kondo_temperature(r, definition="crossing")
    assert kc["T_K"] == kw["T_cross"]
    assert kc["T_K"]/kw["T_K"] == pytest.approx(1.116, abs=1e-3)
    assert kc["chi0"] == pytest.approx(nrgthermo.W_WILSON/(4*kc["T_K"]))


def test_keep_must_be_positive():
    for keep in (0, -1, 2.5, True):
        with pytest.raises(ValueError):
            nrgthermo.kondo_nrg(6, keep=keep)


def test_wilson_definition_warns_with_a_field():
    r = nrgthermo.kondo_nrg(40, J=0.2*A4, Lambda=4.0, keep=200, B=1.2e-6)
    with pytest.warns(UserWarning, match="B = 0"):
        nrgthermo.kondo_temperature(r)


def test_truncated_run_warns():
    """keep = 40 at Lambda = 4 leaves a mean intercept in T chi_imp."""
    r = nrgthermo.kondo_nrg(40, J=0.2*A4, Lambda=4.0, keep=40)
    with pytest.warns(UserWarning, match="raise keep"):
        nrgthermo.kondo_temperature(r)
