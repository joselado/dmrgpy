"""The guard of cvm_solver="variational" (dynamical DMRG).

The variational correction vector minimizes a functional built on the MPO
M = (H-E0-omega)^2 + eta^2, compressed at the DMRG cutoff. On a Kondo
chain with two spins on two Wilson chains (ysr.spin_chain_on_metal) M is
exact on 8 and 12 sites and loses positivity from 16 sites at cutoff
1e-10: <GS|M|GS> is -4.1e-3 against an exact 4.4e-4 on 16 sites, and
-0.152 against 4.0e-4 on 44, where every spectral weight came back
negative while TASK gave 0.0097. v3 and "python" now compare <GS|M|GS>
with ||(H-c)|GS>||^2 + eta^2 from the uncompressed H and raise, and the
CVM driver raises on a negative weight.
"""
import contextlib
import io

import numpy as np
import pytest

from dmrgpy import cppext, cvm, spinchain, ysr

LAMBDA = 2.0
JK = 2*0.15*ysr.lambda_factor(LAMBDA)


def _dimer(n_wilson, version):
    fc = ysr.spin_chain_on_metal(2, n_wilson, site=[0, 1], Js=0.3, JK=JK,
                                 Lambda=LAMBDA, order="sides",
                                 itensor_version=version)
    fc.maxm, fc.nsweeps, fc.cutoff = 64, 10, 1e-10
    fc.cvm_solver, fc.cvm_maxm, fc.cvm_nsweeps = "variational", 64, 2
    return fc


@pytest.mark.skipif(not cppext.available(3), reason="needs the v3 extension")
def test_v3_refuses_a_compressed_operator_that_is_not_positive():
    """16 sites at cutoff 1e-10: M is 11 eta^2 off on the ground state"""
    fc = _dimer(7, 3)
    assert fc.ns == 16
    with contextlib.redirect_stdout(io.StringIO()):
        with pytest.raises(RuntimeError, match="no minimum"):
            ysr.get_stm_spectrum(fc, np.array([0.0]), delta=0.02,
                                 imp=fc.impurities[0], submode="CVM")


@pytest.mark.parametrize("version", [3, "python"])
def test_a_good_operator_passes(version):
    """8 sites: M is exact and the route answers as before"""
    if version == 3 and not cppext.available(3):
        pytest.skip("needs the v3 extension")
    fc = _dimer(3, version)
    assert fc.ns == 8
    with contextlib.redirect_stdout(io.StringIO()):
        y = ysr.get_stm_spectrum(fc, np.array([0.0, 0.3]), delta=0.02,
                                 imp=fc.impurities[0], submode="CVM")[1]
    assert np.all(np.real(y) > 0)


def test_python_refuses_a_compressed_operator_that_is_not_positive():
    """a cutoff far above the default compresses M on 6 spins until it is
    wrong on the ground state; the guard reads the same quantity"""
    sc = spinchain.Spin_Chain(["S=1/2"]*6, itensor_version="python")
    h = sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
            for i in range(5))
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 32, 10
    sc.gs_energy()
    sc.cutoff = 1e-2
    sc.cvm_solver, sc.cvm_maxm, sc.cvm_nsweeps = "variational", 32, 2
    with contextlib.redirect_stdout(io.StringIO()):
        with pytest.raises(RuntimeError, match="no minimum"):
            sc.get_dynamical_correlator(name=(sc.Sz[1], sc.Sz[1]),
                                        es=np.array([0.5]), delta=0.02,
                                        submode="CVM")


def test_a_negative_weight_raises():
    sc = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version="python")
    sc.set_hamiltonian(sum(sc.Sz[i]*sc.Sz[i+1] for i in range(3)))
    sc.gs_energy()
    wf = sc.get_gs()
    cvm._check_ddmrg_values(sc.Sx[0], wf, 0.1, [0.0, 1.0], [0.2, 0.1])
    with pytest.raises(RuntimeError, match="negative spectral weight"):
        cvm._check_ddmrg_values(sc.Sx[0], wf, 0.1, [0.0, 1.0], [0.2, -0.01])
