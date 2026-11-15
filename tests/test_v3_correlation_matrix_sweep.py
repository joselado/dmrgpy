"""The one-sweep correlation matrix of itensor_version=3
(mpscpp3/chain_session.h's Chain::correlation_matrix_sweep, the port of
pyitensor/chain.py's correlation_matrix_sweep) against the per-pair
AutoMPO loop Chain::correlation_matrix, on real, complex, interleaved
spinful and conserved-sector states."""
import warnings

import numpy as np
import pytest

from dmrgpy import cppext, fermionchain, spinchain

warnings.filterwarnings("ignore", message=".*not compiled.*")

pytestmark = pytest.mark.skipif(not cppext.available(3),
                                reason="the itensor_version=3 extension is not compiled")


def _spinless(n, phase=0.0, sector=None):
    fc = fermionchain.Fermionic_Chain(n, itensor_version=3, maxm=40,
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
    # interleaved: two spinless fermion sites per orbital, up then down
    fc = fermionchain.Spinful_Fermionic_Chain(n, itensor_version=3,
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
    "sector_complex": lambda: _spinless(8, phase=0.7, sector=dict(Nf=4)),
}


@pytest.mark.parametrize("case", sorted(CASES))
def test_sweep_is_the_per_pair_loop(case):
    fc = CASES[case]()
    s = fc._session
    ref = s.correlation_matrix(fc.wf0.cpp_handle)
    out = s.correlation_matrix_sweep(fc.wf0.cpp_handle)
    assert out.shape == ref.shape == (fc.ns, fc.ns)
    assert np.max(np.abs(out-ref)) < 1e-12
    assert np.allclose(out, out.conj().T, atol=1e-14)


@pytest.mark.parametrize("nf", [3, 4])
def test_sector_trace_is_the_particle_number(nf):
    fc = _spinless(8, phase=0.7, sector=dict(Nf=nf))
    out = fc._session.correlation_matrix_sweep(fc.wf0.cpp_handle)
    assert np.trace(out).real == pytest.approx(nf, abs=1e-8)


def test_unnormalized_state_scales_as_its_norm_squared():
    fc = _spinless(6, phase=0.4)
    s = fc._session
    wf = fc.wf0*2.0
    ref = s.correlation_matrix(wf.cpp_handle)
    out = s.correlation_matrix_sweep(wf.cpp_handle)
    assert np.max(np.abs(out-ref)) < 1e-12
    assert np.trace(out).real == pytest.approx(4*np.trace(
        s.correlation_matrix_sweep(fc.wf0.cpp_handle)).real, abs=1e-10)


def test_spin_sites_are_refused():
    sc = spinchain.Spin_Chain([2]*4, itensor_version=3, maxm=20, nsweeps=6)
    h = 0
    for i in range(3):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.gs_energy()
    with pytest.raises(ValueError, match="spinless fermion"):
        sc._session.correlation_matrix_sweep(sc.wf0.cpp_handle)


@pytest.mark.parametrize("sector", [None, dict(Nf=3)])
def test_default_takes_the_sweep(sector):
    """get_correlation_matrix() on v3 spinless sites takes dmmode="full",
    now the one sweep, with or without a sector: equal to the per-pair
    AutoMPO loop, and to "fast" outside a sector."""
    from dmrgpy.entropytk import correlationentropy as ce
    fc = _spinless(8, phase=0.7, sector=sector)
    wf = fc.get_gs()
    assert ce._default_dmmode(fc, wf, "electron") == "full"
    got = np.asarray(fc.get_correlation_matrix())
    ref = np.asarray(fc._session.correlation_matrix(wf.cpp_handle))
    assert np.max(np.abs(got-ref)) < 1e-12
    if sector is None:
        fast = np.asarray(fc.get_correlation_matrix(dmmode="fast"))
        assert np.max(np.abs(got-fast)) < 1e-6
