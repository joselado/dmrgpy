"""applyinverse by two-site sweeps (inverse_solver="sweep",
pyitensor/implicit.py's linear_solve, and its port
mpsjulialive/implicit.jl on julia_live), against the exact solve on the
ED matrix of the same chain.

The solver shares its sweep with the implicit imaginary-time factor
(implicit._galerkin_sweep), so test_imaginary_time_implicit.py pins that
the factor did not change.
"""

import warnings

import numpy as np
import pytest

from dmrgpy import cppext, spinchain
from dmrgpy.pyitensor import implicit
from _helpers import julia_live_param

# every julia_live test carries the julia_live id, so -k "not julia_live"
# skips them
JL = pytest.mark.parametrize("version", [julia_live_param()])


def _chain(n=6, maxm=64, solver="sweep", version="python"):
    np.random.seed(1)
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=version)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = maxm
    sc.inverse_solver = solver
    return sc, h


def _exact(sc, A, B):
    """<GS|B^dagger A^-1 B|GS> from the ED matrices"""
    ed = sc.get_ED_obj()
    gs = np.array(sc.get_gs(mode="ED").v)
    b = np.asarray(ed.MO2matrix(B).todense()) @ gs
    Am = np.asarray(ed.MO2matrix(A).todense())
    return np.vdot(b, np.linalg.solve(Am, b))


@pytest.mark.parametrize("which", ["complex", "definite", "indefinite"])
def test_resolvent_matches_ed(which):
    """z - H off the real axis, H - e below the spectrum, and H - e between
    two eigenvalues, at a bond dimension that holds every state"""
    sc, h = _chain()
    e0 = sc.gs_energy()
    E = np.linalg.eigvalsh(np.asarray(sc.get_ED_obj().get_hamiltonian().todense()))
    A = {"complex": (e0 + 1.0 + 0.1j) - h,
         "definite": h - (E[0] - 0.3),
         "indefinite": h - 0.5*(E[3] + E[4])}[which]
    B = sc.Sz[3]
    psi = B*sc.get_gs()
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        x = sc.applyinverse(A, psi, delta=1e-8)
    exact = _exact(sc, A, B)
    assert abs(psi.dot(x) - exact)/abs(exact) < 1e-7


def test_residual_history_falls():
    sc, h = _chain()
    e0 = sc.gs_energy()
    s = sc._session
    A = s._mpo(((e0 + 1.0 + 0.1j) - h).to_terms())
    psi = (sc.Sz[3]*sc.get_gs()).cpp_handle
    x, history = implicit.linear_solve(psi, A, 1e-12, 64, 1e-10, nsweeps=6)
    assert min(history) < 1e-11*history[0]  # below the 1e-8 floor of <r|r>


def test_cvm_explicit_matches_ed():
    sc, h = _chain()
    B = sc.Sz[3]
    es = np.linspace(0.3, 3.0, 6)
    _, d = sc.get_dynamical_correlator(submode="CVM_explicit", name=(B, B),
                                       es=es, delta=0.1)
    _, dex = sc.get_dynamical_correlator(mode="ED", submode="ED", name=(B, B),
                                         es=es, delta=0.1)
    assert np.max(np.abs(np.real(d) - np.real(dex))) < 1e-6*np.max(np.abs(dex))


def test_unknown_solver_raises():
    sc, h = _chain(n=4, solver="gmres")
    psi = sc.Sz[1]*sc.get_gs()
    with pytest.raises(ValueError):
        sc.applyinverse(h - 3.0, psi)


@pytest.mark.skipif(not cppext.available(3),
                    reason="needs the compiled itensor_version=3 extension")
def test_other_backends_raise():
    sc = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version=3)
    sc.set_hamiltonian(sc.Sz[0]*sc.Sz[1] + sc.Sx[1]*sc.Sx[2])
    sc.inverse_solver = "sweep"
    psi = sc.Sz[1]*sc.get_gs()
    with pytest.raises(NotImplementedError):
        sc.applyinverse(sc.Sz[0]*sc.Sz[1] - 3.0, psi)


def test_default_solver_per_backend_and_operator():
    """sweeps for a resolvent off the real axis on "python", BiCGSTAB for a
    Hermitian shift such as the H - e of shift-invert, with a conserved
    sector and on the other backends"""
    from dmrgpy import mpsalgebra
    sc, h = _chain(n=4, solver=None)
    assert sc.inverse_solver is None
    res, shift = (1.0 + 0.1j) - h, h - 0.3
    assert mpsalgebra.resolve_inverse_solver(sc, res) == "sweep"
    assert mpsalgebra.resolve_inverse_solver(sc, shift) == "bicgstab"
    sc.set_conserved_sector(Sz=0)
    assert mpsalgebra.resolve_inverse_solver(sc, res) == "bicgstab"
    sc.set_conserved_sector()
    sc.inverse_solver = "bicgstab"
    assert mpsalgebra.resolve_inverse_solver(sc, res) == "bicgstab"
    sc.inverse_solver = "sweep"
    assert mpsalgebra.resolve_inverse_solver(sc, shift) == "sweep"
    if cppext.available(3):
        sc3 = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version=3)
        assert mpsalgebra.resolve_inverse_solver(sc3, res) == "bicgstab"


@JL
def test_julia_live_resolvent_matches_ed(version):
    """z - H off the real axis by the julia_live sweeps
    (mpsjulialive/implicit.jl), on the chain of test_resolvent_matches_ed"""
    sc, h = _chain(version=version)
    e0 = sc.gs_energy()
    A = (e0 + 1.0 + 0.1j) - h
    B = sc.Sz[3]
    psi = B*sc.get_gs()
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        x = sc.applyinverse(A, psi, delta=1e-8)
    exact = _exact(sc, A, B)
    assert abs(psi.dot(x) - exact)/abs(exact) < 1e-7


@JL
def test_julia_live_default_solver_is_bicgstab(version):
    """inverse_solver=None stays BiCGSTAB on julia_live for every operator,
    as on v3; "sweep" is taken only when asked, and refused in a conserved
    sector"""
    from dmrgpy import mpsalgebra
    sc, h = _chain(n=4, solver=None, version=version)
    res, shift = (1.0 + 0.1j) - h, h - 0.3
    assert mpsalgebra.resolve_inverse_solver(sc, res) == "bicgstab"
    assert mpsalgebra.resolve_inverse_solver(sc, shift) == "bicgstab"
    sc.inverse_solver = "sweep"
    assert mpsalgebra.resolve_inverse_solver(sc, res) == "sweep"
    sc.set_conserved_sector(Sz=0)
    with pytest.raises(ValueError):
        mpsalgebra.resolve_inverse_solver(sc, res)
