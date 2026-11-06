"""cvm_solver="variational" (Jeckelmann's dynamical DMRG) on the compiled
ITensor v3 backend, against the pure-Python backend and against ED.

The chains are small enough that maxm is above the exact bond dimension,
so the only differences left are the solver tolerances. The sweep schedule
is pinned here (cvm_maxm=60, cvm_nsweeps=8) instead of taken from defaults.
"""
import contextlib
import io

import numpy as np
import pytest

from dmrgpy import cppext, cvm, fermionchain, spinchain

pytestmark = pytest.mark.skipif(not cppext.available(3),
                                reason="needs the compiled itensor_version=3 extension")

ES = np.array([-0.5, 0.3, 1.0, 2.2])
ETAS = [0.1, 0.3, 0.5]


def _spin(version, scale=1.0, L=6):
    sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version=version)
    h = 0
    for i in range(L-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(scale*(h + 0.3*sc.Sz[0]))
    return sc, sc.Sz[1], sc.Sz[1]


def _fermion(version, L=6):
    fc = fermionchain.Fermionic_Chain(L, itensor_version=version)
    h = 0
    for i in range(L-1):
        h = h + fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i]
    h = h + 0.4*fc.Cdag[0]*fc.C[0] + 0.5*fc.N[1]*fc.N[2]
    fc.set_hamiltonian(h)
    return fc, fc.C[1], fc.Cdag[1] # A = B^dagger


# measured v3 against python and against ED: 1e-12 to 3e-12 of the peak
# (local CG tolerance 1e-10 on the two-site block); 1e-9 leaves margin
TOL_PY = 1e-9
TOL_ED = 1e-9

MODELS = [("spin", _spin), ("fermion", _fermion)]


def _corr(chain, A, B, eta, mode, submode, es=ES):
    chain.maxm, chain.nsweeps = 60, 12
    chain.cvm_maxm, chain.cvm_nsweeps = 60, 8
    chain.cvm_solver = "variational"
    with contextlib.redirect_stdout(io.StringIO()):
        _x, y = chain.get_dynamical_correlator(name=(A, B), es=es, delta=eta,
                                               mode=mode, submode=submode)
    return np.real(np.asarray(y))


@pytest.mark.parametrize("eta", ETAS)
@pytest.mark.parametrize("name,build", MODELS, ids=[m[0] for m in MODELS])
def test_v3_matches_python_and_ed(name, build, eta):
    chain, A, B = build(3)
    y3 = _corr(chain, A, B, eta, "DMRG", "CVM")
    chain, A, B = build("python")
    yp = _corr(chain, A, B, eta, "DMRG", "CVM")
    chain, A, B = build("python")
    ye = _corr(chain, A, B, eta, "ED", "INV")
    scale = np.max(np.abs(ye))
    assert np.max(np.abs(y3 - yp)) < TOL_PY*scale
    assert np.max(np.abs(y3 - ye)) < TOL_ED*scale


def test_v3_gate_takes_the_variational_route():
    chain, A, B = _spin(3)
    chain.cvm_solver = "variational"
    assert cvm._use_ddmrg(chain, A, B)
    chain.cvm_solver = "cg"
    assert not cvm._use_ddmrg(chain, A, B)


def test_v3_small_scale_hamiltonian():
    """H in units of 1e-3: the solver scale keeps the answer covariant,
    C(eps H; eps omega, eps eta) = C(H; omega, eta)/eps."""
    eps = 1e-3
    chain, A, B = _spin("python")
    ref = _corr(chain, A, B, 0.3, "ED", "INV")
    chain, A, B = _spin(3, scale=eps)
    y = _corr(chain, A, B, 0.3*eps, "DMRG", "CVM", es=ES*eps)
    assert np.max(np.abs(y*eps - ref)) < TOL_ED*np.max(np.abs(ref))


def test_v3_two_sites_route_to_ed():
    """v3 with fewer than 3 sites goes to ED (mode.py) before the sweep."""
    sc = spinchain.Spin_Chain(["S=1/2"]*2, itensor_version=3)
    sc.set_hamiltonian(sc.Sx[0]*sc.Sx[1] + sc.Sy[0]*sc.Sy[1] + sc.Sz[0]*sc.Sz[1])
    y = _corr(sc, sc.Sz[0], sc.Sz[0], 0.3, "DMRG", "CVM")
    ref = _corr(sc, sc.Sz[0], sc.Sz[0], 0.3, "ED", "INV")
    assert np.max(np.abs(y - ref)) < TOL_ED*np.max(np.abs(ref))
