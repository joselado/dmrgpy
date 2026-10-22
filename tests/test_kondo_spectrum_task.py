"""get_kondo_spectrum(mode="DMRG", submode="TASK"): the exact sums of the
ED route on the tangent-space Krylov spectrum (kondospectrumtk/taskkondo.py).

With cutoff=0 and a bond dimension that saturates the Hilbert space, the
tangent space is the whole complement of the ground state, so once the
block Krylov space is exhausted the spectrum must be the ED one at T=0,
second order, third-order Kondo and potential term alike.
"""

import numpy as np
import pytest

from dmrgpy import cppext, spinchain

VERSIONS = ["python", pytest.param(3, marks=pytest.mark.skipif(
    not cppext.available(3), reason="v3 extension not compiled"))]

G, MUB = 2.0, 5.7883818066e-5
EV = np.linspace(-6e-3, 6e-3, 61)
KW = dict(site=0, Jrho_s=-0.1, U=0.3, T=0., omega0=20e-3, Gamma0=5e-5, order=3)


def impurity_chain(nb=4, maxm=64, version="python"):
    """An S=1 impurity with easy-plane anisotropy in a field, coupled to a
    short S=1/2 chain."""
    sc = spinchain.Spin_Chain(["S=1"] + ["S=1/2"] * nb, itensor_version=version)
    h = 1e-3 * sc.Sz[0] * sc.Sz[0] + G * MUB * 2.0 * sc.Sz[0]
    for i in range(nb):
        J = 0.6e-3 if i == 0 else 1e-3
        h = h + J * (sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1]
                     + sc.Sz[i] * sc.Sz[i + 1])
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps, sc.cutoff = maxm, 20, 0.
    return sc


@pytest.mark.parametrize("version", VERSIONS)
def test_saturated_task_matches_ed_at_zero_temperature(version):
    sc = impurity_chain(version=version)
    _, ed = sc.get_kondo_spectrum(EV, mode="ED", **KW)
    _, tk = sc.get_kondo_spectrum(EV, mode="DMRG", submode="TASK", nkry=50, **KW)
    assert np.max(np.abs(tk - ed)) < 1e-9 * np.max(np.abs(ed))


def test_task_rejects_the_grid_keywords_of_the_other_routes():
    sc = impurity_chain()
    with pytest.raises(TypeError, match="es"):
        sc.get_kondo_spectrum(EV, mode="DMRG", submode="TASK",
                              es=np.linspace(-1e-2, 1e-2, 10), **KW)
    with pytest.raises(TypeError, match="delta"):
        sc.get_kondo_spectrum(EV, mode="DMRG", submode="TASK", delta=1e-5, **KW)


@pytest.mark.parametrize("version", VERSIONS)
def test_task_is_the_default(version):
    sc = impurity_chain(version=version)
    _, ed = sc.get_kondo_spectrum(EV, mode="ED", **KW)
    _, tk = sc.get_kondo_spectrum(EV, mode="DMRG", nkry=50, **KW)
    assert np.max(np.abs(tk - ed)) < 1e-9 * np.max(np.abs(ed))


@pytest.mark.parametrize("version", VERSIONS)
def test_task_on_a_chain_too_short_for_the_session_is_the_ed_spectrum(version):
    """v3 below 3 sites and "python" on 1 site run on ED (mode.py), and so
    does the TASK route, as every other DMRG route does."""
    sc = spinchain.Spin_Chain(["S=1/2"], itensor_version=version)
    sc.set_hamiltonian(G * MUB * 5.0 * sc.Sz[0])
    kw = dict(KW, U=0.1)
    _, ed = sc.get_kondo_spectrum(EV, mode="ED", **kw)
    _, tk = sc.get_kondo_spectrum(EV, mode="DMRG", **kw)
    assert np.max(np.abs(tk - ed)) < 1e-12 * np.max(np.abs(ed))


def test_task_needs_a_backend_that_has_it():
    sc = spinchain.Spin_Chain(["S=1/2"] * 3, itensor_version=2)
    sc.set_hamiltonian(sc.Sx[0] * sc.Sx[1] + sc.Sz[1] * sc.Sz[2])
    if sc.get_mode(mode="DMRG") != "DMRG":
        pytest.skip("extension not compiled")
    with pytest.raises(NotImplementedError, match="python"):
        sc.get_kondo_spectrum(EV, mode="DMRG", submode="TASK", **KW)


@pytest.mark.parametrize("version", VERSIONS)
def test_task_averages_a_kramers_doublet_like_ed(version):
    """A U(1)-symmetric 3-site chain has a ground doublet related only by
    spin flip, so a single member is not the average, n_gs=2 is. The
    partner state comes out of the Ritz solve at +-1e-16 and has to be
    snapped onto the zero-bias step, or eV=0 missed ED by 37% in 5 runs of
    6 (TaskKondoSpectrum's degeneracy_tol). Repeated, since the members
    DMRG returns differ run to run."""
    sc = spinchain.Spin_Chain(["S=1/2"] * 3, itensor_version=version)
    h = 0
    for i in range(2):
        h = h + 1e-3 * (sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1]
                        + sc.Sz[i] * sc.Sz[i + 1])
    h = h + 0.3e-3 * sc.Sz[0] * sc.Sz[2]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps, sc.cutoff = 8, 20, 0.
    eV = np.linspace(-4e-3, 4e-3, 81)  # eV=0 is on the grid
    kw = dict(site=0, Jrho_s=-0.1, U=0., T=0., omega0=20e-3, Gamma0=5e-5, order=3)
    _, ed = sc.get_kondo_spectrum(eV, mode="ED", **kw)
    for _ in range(3):
        sc.restart()
        _, tk = sc.get_kondo_spectrum(eV, mode="DMRG", submode="TASK", nkry=30,
                                      n_gs=2, **kw)
        assert np.max(np.abs(tk - ed)) < 1e-9 * np.max(np.abs(ed))
