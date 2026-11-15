"""The one-sweep correlation matrix of itensor_version="python"
(pyitensor/chain.py's correlation_matrix_sweep, the algorithm of
ITensorMPS.jl's correlation_matrix) against the per-pair MPO loop it
replaced, on real, complex, spinful and sector states."""
import numpy as np
import pytest

from dmrgpy import fermionchain


def _spinless(n, phase=0.0, sector=None):
    fc = fermionchain.Fermionic_Chain(n, itensor_version="python", maxm=40,
                                      nsweeps=8)
    h = 0
    for i in range(n-1):
        h = h + np.exp(1j*phase)*fc.Cdag[i]*fc.C[i+1]
    h = h + h.get_dagger()
    h = h + sum(1.5*fc.N[i]*fc.N[i+1] for i in range(n-1)) + 0.2*fc.N[0]
    fc.set_hamiltonian(h)
    if sector: fc.set_conserved_sector(**sector)
    fc.gs_energy()
    return fc


def _spinful(n, phase=0.3):
    fc = fermionchain.Spinful_Fermionic_Chain(n, itensor_version="python",
                                              maxm=40, nsweeps=8)
    h = 0
    for i in range(n-1):
        h = h + np.exp(1j*phase)*(fc.Cdagup[i]*fc.Cup[i+1]
                                  + fc.Cdagdn[i]*fc.Cdn[i+1])
    h = h + h.get_dagger()
    h = h + 2.0*sum(fc.Nup[i]*fc.Ndn[i] for i in range(n)) + 0.3*fc.Nup[0]
    fc.set_hamiltonian(h)
    fc.gs_energy()
    return fc


CASES = {
    "real": lambda: _spinless(8),
    "complex": lambda: _spinless(8, phase=0.7),
    "spinful": lambda: _spinful(4),
    "sector": lambda: _spinless(8, sector=dict(Nf=3)),
}


@pytest.mark.parametrize("case", sorted(CASES))
def test_sweep_is_the_per_pair_loop(case):
    fc = CASES[case]()
    s = fc._session
    ref = s.correlation_matrix_per_pair(fc.wf0.cpp_handle)
    out = s.correlation_matrix(fc.wf0.cpp_handle)
    assert np.max(np.abs(out-ref)) < 1e-12
    assert np.allclose(out, out.conj().T, atol=1e-14)


def test_unnormalized_state_scales_as_its_norm_squared():
    fc = _spinless(6)
    s = fc._session
    wf = fc.wf0*2.0
    ref = s.correlation_matrix_per_pair(wf.cpp_handle)
    out = s.correlation_matrix(wf.cpp_handle)
    assert np.max(np.abs(out-ref)) < 1e-12
    assert np.trace(out).real == pytest.approx(4*np.trace(
        s.correlation_matrix(fc.wf0.cpp_handle)).real, abs=1e-10)


def test_sweep_runs_on_jax():
    jax = pytest.importorskip("jax")
    from dmrgpy.pyitensor import backend
    fc = _spinless(6, phase=0.4)
    s = fc._session
    ref = s.correlation_matrix_per_pair(fc.wf0.cpp_handle)
    try:
        backend.set_backend("jax")
        out = s.correlation_matrix(fc.wf0.cpp_handle)
    finally:
        backend.set_backend("numpy")
    assert np.max(np.abs(np.asarray(out)-ref)) < 1e-10
