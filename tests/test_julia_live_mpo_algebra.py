"""julia_live MPOs: parity-odd fermion operators (a single C or Cdag).

ITensorMPS's OpSum refuses a parity-odd fermionic term, so
mpsjulialive/mpo.py sends such an operator in the Jordan-Wigner-dressed
form MultiOperator.to_terms() gives (A/Adag/F names, bosonic for OpSum,
see mpsjulialive/fermion_ops.jl). Everything is checked against
mode="ED".
"""
import numpy as np
import pytest

from dmrgpy import fermionchain, spinchain, timedependent
from _helpers import julia_live_param, setup_backend

BACKENDS = [julia_live_param(), pytest.param("python", id="python")]


def _fermion_chain(itensor_version, n=4):
    """Spinless fermions with first-neighbor hopping and a weak staggered
    onsite energy, so the ground state is not degenerate"""
    fc = fermionchain.Fermionic_Chain(n)
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i]
    for i in range(n):
        h = h + 0.2*(-1)**i*fc.N[i] + 0.1*i*fc.N[i]
    fc.set_hamiltonian(h)
    setup_backend(fc, itensor_version)
    fc.maxm = 30
    fc.nsweeps = 10
    return fc


def _sp(sc, i):
    """S+ on site i"""
    return sc.Sx[i] + 1j*sc.Sy[i]


def _spin_chain(itensor_version, n=4):
    """Heisenberg chain with a field on the first two sites only, which
    breaks the SU(2) symmetry, so that <S+> is not zero"""
    sc = spinchain.Spin_Chain(["1/2"]*n)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    h = h + 0.6*sc.Sx[0] + 0.3*sc.Sz[1]
    sc.set_hamiltonian(h)
    setup_backend(sc, itensor_version)
    sc.maxm = 30
    sc.nsweeps = 10
    return sc


@pytest.mark.parametrize("itensor_version", BACKENDS)
def test_single_fermion_operators_reproduce_ed(itensor_version):
    """<Cdag_i C_j> three ways: the even MultiOperator, the product of the
    two odd MPOs, and the overlap of two odd operators applied to the
    ground state. The pair (2,0) checks the sign of the string."""
    fc = _fermion_chain(itensor_version)
    gs = fc.get_gs()
    for (i, j) in [(0, 2), (2, 0), (1, 2), (0, 3)]:
        ref = fc.vev(fc.Cdag[i]*fc.C[j], mode="ED")
        assert abs(ref) > 1e-2 # a pair that is not zero by symmetry
        assert fc.vev(fc.Cdag[i]*fc.C[j]) == pytest.approx(ref, abs=1e-6)
        ov = (fc.C[i]*gs).dot(fc.C[j]*gs)
        assert ov == pytest.approx(ref, abs=1e-6)


@pytest.mark.parametrize("submode", ["KPM", "CVM"])
def test_fermion_dynamical_correlator_julia_matches_ed(submode):
    """The single-particle Green's function on julia_live, which needs the
    odd MPOs of Cdag_0 and C_2, against mode="ED". CVM agrees with ED to
    1e-6. KPM on every DMRG backend differs from the ED KPM by about 2e-4
    here, from its own estimate of the spectral window, so julia_live KPM
    is pinned to "python" KPM to 1e-6 and to ED to 1e-3 (a wrong string
    sign would move it by the height of the peak, 0.5)."""
    pytest.importorskip("juliacall")
    fc = _fermion_chain("julia_live")
    fp = _fermion_chain("python")
    es = np.linspace(-3., 3., 13)
    for (A, B) in [(0, 2), (1, 1)]:
        for (Aname, Bname) in [("Cdag", "C"), ("C", "Cdag")]:
            kw = dict(es=es, delta=0.2, submode=submode)
            name = (getattr(fc, Aname)[A], getattr(fc, Bname)[B])
            x, y = fc.get_dynamical_correlator(mode="DMRG", name=name, **kw)
            xe, ye = fc.get_dynamical_correlator(mode="ED", name=name, **kw)
            assert np.allclose(x, xe)
            assert np.max(np.abs(ye)) > 1e-2
            tol = 1e-6 if submode == "CVM" else 1e-3
            assert np.max(np.abs(np.array(y) - np.array(ye))) < tol
            if submode == "KPM":
                namep = (getattr(fp, Aname)[A], getattr(fp, Bname)[B])
                _xp, yp = fp.get_dynamical_correlator(mode="DMRG",
                                                      name=namep, **kw)
                assert np.max(np.abs(np.array(y) - np.array(yp))) < 1e-6


def test_fermion_evolution_aba_julia_matches_ed():
    """<GS|C_0 e^{iHt} N_2 e^{-iHt} Cdag_0|GS> on julia_live against ED"""
    pytest.importorskip("juliacall")
    fc = _fermion_chain("julia_live")
    kw = dict(nt=40, dt=0.05, A=fc.Cdag[0], B=fc.N[2])
    _t, y = timedependent.evolution_ABA(fc, mode="DMRG", **kw)
    _te, ye = timedependent.evolution_ABA(fc, mode="ED", **kw)
    assert np.max(np.abs(np.array(y) - np.array(ye))) < 1e-5
