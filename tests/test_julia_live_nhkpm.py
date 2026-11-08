"""NH-KPM (nonhermitian/kpm.py) on itensor_version="julia_live": the
coupled forward/adjoint Chebyshev recursion of mpsjulialive/nhkpm.jl,
against the ED route and "python" at the same E_max and n. The chain is
the asymmetric-hopping one of test_nh_dmrg.py, non-symmetric as well as
non-Hermitian, so a left vector taken as the transpose instead of the
adjoint would show, with two onsite terms that make its biorthogonal
ground state unique. NH-DMRG starts from a random MPS, so the julia_live
solve is repeated on two fresh chains."""
import numpy as np
import pytest

from dmrgpy import fermionchain
from dmrgpy.nonhermitian import kpm as nhkpm

from _helpers import julia_live_param, setup_backend

KW = dict(es=np.linspace(0.0, 4.0, 5), delta=0.3, E_max=10.0, n=40)


def chain(version, n=4, tr=1.0, tl=0.6, gamma=0.7):
    fc = fermionchain.Fermionic_Chain(n)
    setup_backend(fc, version)
    h = 0
    for i in range(n - 1):
        h = h + tr * fc.Cdag[i] * fc.C[i + 1] + tl * fc.Cdag[i + 1] * fc.C[i]
    for i in range(n):
        h = h + 1j * gamma * (-1)**i * fc.Cdag[i] * fc.C[i]
    for i in range(n - 1):
        h = h + 0.5 * (fc.N[i] - 0.5) * (fc.N[i + 1] - 0.5)
    # breaks the PT symmetry that ties a conjugate pair at the smallest
    # real part, so that the biorthogonal ground state is unique and the
    # ED route can pair its left and right vectors
    h = h + 0.25j * fc.N[0] + 0.15 * fc.N[2]
    fc.set_hamiltonian(h)
    fc.maxm, fc.nsweeps = 20, 20
    fc.kpmmaxm, fc.kpmcutoff = 64, 1e-14
    return fc


def pairs(fc):
    """A diagonal pair, and a pair of hoppings across the chain, whose
    Jordan-Wigner strings cross two sites, where a wrong sign would show.
    Both are parity-even: julia_live builds no MPO of a single C yet."""
    return [(fc.N[1], fc.N[1]),
            (fc.Cdag[0] * fc.C[3], fc.Cdag[3] * fc.C[1])]


def ed_reference(pair_index):
    fc = chain("python")
    A, B = pairs(fc)[pair_index]
    _, y = fc.get_dynamical_correlator(mode="ED", submode="KPM",
                                       name=(A, B), **KW)
    return np.asarray(y)


@pytest.mark.parametrize("version", ["python", julia_live_param()])
@pytest.mark.parametrize("pair_index", [0, 1])
def test_nhkpm_matches_ed(version, pair_index):
    ref = ed_reference(pair_index)
    assert np.max(np.abs(ref)) > 1e-3
    for _ in range(2 if version == "julia_live" else 1):
        fc = chain(version)
        A, B = pairs(fc)[pair_index]
        _, y = nhkpm.dynamical_correlator_nhkpm(fc, name=(A, B), **KW)
        assert np.max(np.abs(np.asarray(y) - ref)) < 1e-6 * np.max(np.abs(ref))


@pytest.mark.parametrize("version", [julia_live_param()])
def test_nhkpm_julia_matches_python(version):
    py = chain("python")
    A, B = pairs(py)[1]
    _, yp = nhkpm.dynamical_correlator_nhkpm(py, name=(A, B), **KW)
    jl = chain(version)
    A, B = pairs(jl)[1]
    _, yj = nhkpm.dynamical_correlator_nhkpm(jl, name=(A, B), **KW)
    assert np.max(np.abs(np.asarray(yj) - yp)) < 1e-6 * np.max(np.abs(yp))


@pytest.mark.parametrize("version", [julia_live_param()])
def test_nhkpm_julia_refuses_a_set_state(version):
    fc = chain(version)
    fc.gs_energy()
    fc.set_gs(fc.get_gs().copy())
    with pytest.raises(RuntimeError, match="left eigenvector"):
        nhkpm.dynamical_correlator_nhkpm(fc, name=pairs(fc)[0], **KW)
