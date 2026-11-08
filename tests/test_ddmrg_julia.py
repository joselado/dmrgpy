"""cvm_solver="variational" (Jeckelmann's dynamical DMRG) on
itensor_version="julia_live" (mpsjulialive/ddmrg.jl), against the
pure-Python backend and against ED.

The chains are small enough that cvm_maxm is above the exact bond
dimension, so the only differences left are the solver tolerances. The
sweep schedule is pinned here (cvm_maxm=60, cvm_nsweeps=8), as in
test_ddmrg_v3.py.
"""
import contextlib
import io

import numpy as np
import pytest

from dmrgpy import cppext, cvm, spinchain
from _helpers import julia_live_param

# every test carries the julia_live id, so -k "not julia_live" skips them
JL = pytest.mark.parametrize("version", [julia_live_param()])

ES = np.array([-0.5, 0.3, 1.0, 2.2])
ETAS = [0.1, 0.3]

# measured julia_live against python and against ED on these chains:
# 2e-14 to 2e-13 of the peak (local CG tolerance 1e-10 on the two-site
# block); 1e-9 leaves margin
TOL = 1e-9


def _spin(version, scale=1.0, L=6):
    sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version=version)
    h = 0
    for i in range(L-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(scale*(h + 0.3*sc.Sz[0]))
    return sc


def _corr(chain, A, B, eta, mode, submode, es=ES, solver="variational"):
    chain.maxm, chain.nsweeps = 60, 12
    chain.cvm_maxm, chain.cvm_nsweeps = 60, 8
    chain.cvm_solver = solver
    with contextlib.redirect_stdout(io.StringIO()):
        _x, y = chain.get_dynamical_correlator(name=(A, B), es=es, delta=eta,
                                               mode=mode, submode=submode)
    return np.real(np.asarray(y))


@pytest.fixture
def julia_calls(monkeypatch):
    """Count the calls into ddmrg.jl, so a test can tell the variational
    Julia route from the CG one (no backend substitution)."""
    from dmrgpy.mpsjulialive import ddmrg as ddmrgjl
    calls = []
    main = ddmrgjl.Mainjl # the real Julia Main, captured before the patch
    real = main.ddmrg_correction_vector

    class _Spy:
        def __getattr__(self, name):
            return getattr(main, name)

        def ddmrg_correction_vector(self, *args):
            calls.append(args[3])
            return real(*args)
    monkeypatch.setattr(ddmrgjl, "Mainjl", _Spy())
    return calls


@JL
@pytest.mark.parametrize("L", [6, 8])
def test_julia_matches_python_and_ed(julia_calls, version, L):
    for eta in ETAS:
        chain = _spin(version, L=L)
        yj = _corr(chain, chain.Sz[1], chain.Sz[1], eta, "DMRG", "CVM")
        chain = _spin("python", L=L)
        yp = _corr(chain, chain.Sz[1], chain.Sz[1], eta, "DMRG", "CVM")
        ye = _corr(chain, chain.Sz[1], chain.Sz[1], eta, "ED", "INV")
        scale = np.max(np.abs(ye))
        assert np.max(np.abs(yj - yp)) < TOL*scale, eta
        assert np.max(np.abs(yj - ye)) < TOL*scale, eta
        if cppext.available(3): # and v3's own sweep, where it is compiled
            chain = _spin(3, L=L)
            y3 = _corr(chain, chain.Sz[1], chain.Sz[1], eta, "DMRG", "CVM")
            assert np.max(np.abs(yj - y3)) < TOL*scale, eta
    assert len(julia_calls) == len(ES)*len(ETAS) # every point ran Julia


@JL
def test_julia_small_scale_hamiltonian(julia_calls, version):
    """H in units of 1e-3: the solver scale keeps the answer covariant,
    C(eps H; eps omega, eps eta) = C(H; omega, eta)/eps."""
    eps = 1e-3
    chain = _spin("python")
    ref = _corr(chain, chain.Sz[1], chain.Sz[1], 0.3, "ED", "INV")
    chain = _spin(version, scale=eps)
    y = _corr(chain, chain.Sz[1], chain.Sz[1], 0.3*eps, "DMRG", "CVM",
              es=ES*eps)
    assert len(julia_calls) == len(ES)
    assert np.max(np.abs(y*eps - ref)) < TOL*np.max(np.abs(ref))


@JL
def test_julia_non_adjoint_pair_goes_to_cg(julia_calls, version):
    """(Sz0, Sz3) has no variational functional: the gate sends it to the
    CG, which still matches ED, and ddmrg.jl is never called."""
    chain = _spin(version)
    chain.cvm_solver = "variational"
    A, B = chain.Sz[0], chain.Sz[3]
    assert cvm._use_ddmrg(chain, A, A)
    assert not cvm._use_ddmrg(chain, A, B)
    es = ES[1:3]
    y = _corr(chain, A, B, 0.3, "DMRG", "CVM", es=es)
    ref = _corr(chain, A, B, 0.3, "ED", "INV", es=es)
    assert julia_calls == []
    assert np.max(np.abs(y - ref)) < 1e-4*np.max(np.abs(ref))
    chain.cvm_solver = "cg"
    assert not cvm._use_ddmrg(chain, A, A)


@JL
def test_julia_keeps_ground_state_when_cvm_maxm_differs(julia_calls, version):
    """maxm != cvm_maxm: the variational route must use the chain's own
    ground state and leave it in place, not re-solve it at cvm_maxm."""
    chain = _spin(version, L=8)
    chain.maxm, chain.nsweeps = 40, 12
    chain.gs_energy()
    wf0, e0 = chain.get_gs(), chain.e0
    chain.cvm_maxm, chain.cvm_nsweeps = 60, 8
    chain.cvm_solver = "variational"
    with contextlib.redirect_stdout(io.StringIO()):
        _x, y = chain.get_dynamical_correlator(name=(chain.Sz[1], chain.Sz[1]),
                    es=ES, delta=0.3, mode="DMRG", submode="CVM")
    y = np.real(np.asarray(y))
    assert len(julia_calls) == len(ES)
    assert chain.wf0 is wf0
    assert chain.e0 == e0
    assert chain.maxm == 40
    assert chain.get_gs() is wf0 # still current at the restored maxm
    pc = _spin("python", L=8)
    ref = _corr(pc, pc.Sz[1], pc.Sz[1], 0.3, "ED", "INV")
    assert np.max(np.abs(y - ref)) < TOL*np.max(np.abs(ref))
