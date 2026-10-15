"""Static observables of infinite chains (2026-10-15).

Two findings from the user-guide figure runs.

A. iDMRG on the S=1/2 Heisenberg chain gave alternating bonds: on n_uc=1 at
maxm=16, 3<Sz0 Sz1> = -0.4137 against e = -0.4431; on n_uc=2 the two bonds
were -0.4719 and -0.4143. The suspected missing lambda^-1 in the cell
extraction is not there: the growing chain's own central bond alternates
between -0.4721 and -0.4140 from one micro-step to the next, the extracted
cell reproduces both, the mean bond equals e on n_uc=2 to the truncation
error, and models whose cell is exact (AKLT at maxm=2, Majumdar-Ghosh on a
two-site cell) give exact tiled correlators. The iDMRG fixed point at
finite maxm is dimerized, and the alternation shrinks as maxm grows. On
n_uc=1 a one-site observable reads only the first bond of the two-site
cell, so gs_energy, vev and correlator now warn, with the entanglement
entropies of the two bonds.

B. VUMPS states with converged=False gave correlators about 10 per cent
off while the energy looked right. vev, correlator and the excitation
ansatz now warn with the gauge mismatch.
"""
import warnings

import numpy as np
import pytest

from dmrgpy import cppext
from dmrgpy import infinitechain
from dmrgpy.pyitensor import idmrg

BACKENDS = ["python",
            pytest.param(3, marks=pytest.mark.skipif(
                not cppext.available(3), reason="needs the mpscpp3 extension"))]
OPS = ("Sx", "Sy", "Sz")


def _dot(ic, ga, i, gb, j):
    return sum(getattr(ic, o + ga)[i]*getattr(ic, o + gb)[j] for o in OPS)


def _bond(ic, p, r=1):
    return sum(ic.correlator(o, p, o, r).real for o in OPS)


def _idmrg(itv, sites, h_of, maxm, maxiter=200):
    ic = infinitechain.Infinite_Spin_Chain(sites, itensor_version=itv)
    ic.set_hamiltonian(h_of(ic))
    ic.gs_method = "idmrg"
    ic.maxm, ic.maxiter = maxm, maxiter
    return ic


def _heisenberg(n_uc):
    def h_of(ic):
        if n_uc == 1:
            return _dot(ic, "C", 0, "R", 0)
        return _dot(ic, "C", 0, "C", 1) + _dot(ic, "C", 1, "R", 0)
    return h_of


def _no_warning(fn):
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        return fn()


@pytest.mark.parametrize("itv", BACKENDS)
def test_aklt_tiled_cell_is_exact_and_uniform(itv):
    def h_of(ic):
        d = _dot(ic, "C", 0, "R", 0)
        return d + d*d*(1/3.)
    ic = _idmrg(itv, ["1"], h_of, maxm=2, maxiter=100)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        e = ic.gs_energy()
        c = [ic.correlator("Sz", 0, "Sz", r).real for r in (1, 2, 3)]
    assert not [w for w in caught if "translation invariant" in str(w.message)]
    assert e == pytest.approx(-2/3, abs=1e-9)
    for r in (1, 2, 3):
        assert c[r - 1] == pytest.approx((4/3)*(-1/3)**r, abs=1e-7)
    assert ic._cell_asymmetry[0] < 1e-4


@pytest.mark.parametrize("itv", BACKENDS)
def test_majumdar_ghosh_two_site_cell_is_exact(itv):
    def h_of(ic):
        return (_dot(ic, "C", 0, "C", 1) + _dot(ic, "C", 1, "R", 0)
                + 0.5*(_dot(ic, "C", 0, "R", 0) + _dot(ic, "C", 1, "R", 1)))
    ic = _idmrg(itv, ["1/2", "1/2"], h_of, maxm=2, maxiter=100)
    e = ic.gs_energy()
    assert e == pytest.approx(-3/8, abs=1e-12)
    b = sorted([_bond(ic, 0), _bond(ic, 1)])
    assert b[0] == pytest.approx(-0.75, abs=1e-10)
    assert b[1] == pytest.approx(0.0, abs=1e-10)
    nnn = [_bond(ic, p, 2) for p in (0, 1)]
    assert np.allclose(nnn, 0.0, atol=1e-10)


@pytest.mark.parametrize("itv", BACKENDS)
def test_heisenberg_two_site_cell_mean_bond_is_the_energy(itv):
    """The bonds alternate (finite-maxm dimerization of the iDMRG fixed
    point), and their mean is e to the truncation error: 4.2e-7 at maxm=16
    on both backends, unchanged when etol goes from 1e-10 to 1e-13."""
    ic = _idmrg(itv, ["1/2", "1/2"], _heisenberg(2), maxm=16)
    e = ic.gs_energy()
    b0, b1 = _bond(ic, 0), _bond(ic, 1)
    assert 0.5*(b0 + b1) == pytest.approx(e, abs=2e-6)
    assert abs(b0 - b1) > 0.03
    assert e == pytest.approx(0.25 - np.log(2), abs=1e-4)


@pytest.mark.parametrize("itv", BACKENDS)
def test_heisenberg_one_site_cell_warns_that_it_is_dimerized(itv):
    ic = _idmrg(itv, ["1/2"], _heisenberg(1), maxm=16)
    with pytest.warns(RuntimeWarning, match="not translation invariant"):
        e = ic.gs_energy()
    dist, s0, s1 = ic._cell_asymmetry
    assert dist > 0.5 and abs(s0 - s1) > 0.03
    with pytest.warns(RuntimeWarning, match="read only the first bond"):
        b0 = _bond(ic, 0)
    with pytest.warns(RuntimeWarning, match="not translation invariant"):
        ic.vev("Sz", 0)
    # the first bond of the cell is the weak one, far from e
    assert b0 - e > 0.02


def test_one_site_cell_both_bonds_average_to_the_energy():
    """On "python" the second bond of the n_uc=1 cell can be read directly:
    the two bonds are the n_uc=2 pair, and their mean is e (the residual
    1.3e-5 is the unconverged alternation of the density itself)."""
    ic = _idmrg("python", ["1/2"], _heisenberg(1), maxm=16)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        e = ic.gs_energy()
    b = [sum(idmrg.two_point_correlator(ic._result, o, p, o, 1).real
             for o in OPS) for p in (0, 1)]
    assert np.mean(b) == pytest.approx(e, abs=1e-4)
    assert abs(b[0] - b[1]) > 0.03


@pytest.mark.skipif(not cppext.available(3), reason="needs the mpscpp3 extension")
def test_cell_bond_weights_agree_between_backends():
    w = []
    for itv in ("python", 3):
        ic = _idmrg(itv, ["1/2"], _heisenberg(1), maxm=16)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            ic.gs_energy()
        w.append(ic._cell_asymmetry)
    assert np.allclose(w[0], w[1], atol=1e-5)


def _unconverged_vumps(itv):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=itv)
    ic.set_hamiltonian(_heisenberg(1)(ic))
    ic.gs_method = "vumps"
    ic.maxm, ic.maxiter, ic.vumps_nrestarts = 4, 20, 1
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        ic.gs_energy()
    return ic


@pytest.mark.parametrize("itv", BACKENDS)
def test_unconverged_vumps_observables_warn_with_the_mismatch(itv):
    ic = _unconverged_vumps(itv)
    assert ic.converged is False
    with pytest.warns(RuntimeWarning, match="did not converge.*gauge mismatch"):
        ic.vev("Sz", 0)
    with pytest.warns(RuntimeWarning, match="did not converge.*gauge mismatch"):
        ic.correlator("Sz", 0, "Sz", 1)


@pytest.mark.parametrize("itv", BACKENDS)
def test_converged_vumps_observables_do_not_warn(itv):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=itv)
    ic.set_hamiltonian(-4*ic.SzC[0]*ic.SzR[0] - 3.0*ic.SxC[0])
    ic.gs_method = "vumps"
    ic.maxm, ic.maxiter = 4, 300
    ic.gs_energy()
    assert ic.converged
    _no_warning(lambda: ic.vev("Sx", 0))
    _no_warning(lambda: ic.correlator("Sz", 0, "Sz", 1))
