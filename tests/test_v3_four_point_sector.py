"""The four-point tensor <Cdag_i C_j Cdag_k C_l> on itensor_version=3 in a
conserved sector. The one-sweep and fold routes closed their right end
with a link of the wrong QN arrow, which ITensor aborts the process on
("Mismatched QN Index arrows"); the per-tuple AutoMPO routes
(ctmode="full") abort inside ITensor on some tuples and now raise
ValueError instead."""
import warnings

import numpy as np
import pytest

from dmrgpy import cppext, fermionchain

pytestmark = pytest.mark.skipif(not cppext.available(3),
                                reason="no compiled mpscpp3 extension")


def _spinless(sector=None, n=6):
    fc = fermionchain.Fermionic_Chain(n, itensor_version=3, maxm=40,
                                      nsweeps=10)
    h = sum(0.7*np.exp(0.4j)*fc.Cdag[i]*fc.C[i+1] for i in range(n-1))
    h = h + h.get_dagger() + sum(fc.N[i]*fc.N[i+1] for i in range(n-1))
    fc.set_hamiltonian(h)
    if sector:
        fc.set_conserved_sector(**sector)
    return fc


def _native(sector=None, n=4):
    fc = fermionchain.Spinful_Fermionic_Chain_Native(n, itensor_version=3)
    h = 0
    for i in range(n-1):
        h = h + (0.9+0.2j)*fc.Cdagup[i]*fc.Cup[i+1]
        h = h + (0.9+0.2j)*fc.Cdagdn[i]*fc.Cdn[i+1]
    for i in range(n):
        h = h + 1.3*(fc.Nup[i]-.5)*(fc.Ndn[i]-.5)
    h = h + h.get_dagger() + 0.3*(fc.Nup[0]-fc.Ndn[0])
    fc.maxm = 30
    fc.nsweeps = 12
    fc.set_hamiltonian(h)
    if sector:
        fc.set_conserved_sector(**sector)
    return fc


def _tensor(wf, ctmode=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if ctmode is None:
            return np.asarray(wf.get_four_correlation_tensor())
        return np.asarray(wf.get_four_correlation_tensor(ctmode=ctmode))


def test_spinless_sweep_in_a_sector_is_the_explicit_loop():
    wf = _spinless(dict(Nf=3)).get_gs()
    ref = _tensor(wf, "explicit")
    assert np.max(np.abs(ref)) > 1e-2
    assert np.max(np.abs(_tensor(wf)-ref)) < 1e-10
    assert np.max(np.abs(_tensor(wf, "sweep")-ref)) < 1e-10


def test_spinless_full_in_a_sector_raises():
    wf = _spinless(dict(Nf=3)).get_gs()
    with pytest.raises(ValueError, match="sweep"):
        _tensor(wf, "full")


def test_spinless_without_a_sector_is_unchanged():
    wf = _spinless().get_gs()
    ref = _tensor(wf, "explicit")
    for ctmode in (None, "sweep", "full"):
        assert np.max(np.abs(_tensor(wf, ctmode)-ref)) < 1e-10


def test_native_spinful_in_an_nf_sector():
    wf = _native(dict(Nf=4)).get_gs()
    ref = _tensor(wf, "explicit")
    assert np.max(np.abs(ref)) > 1e-2
    for ctmode in (None, "fold", "full"):
        assert np.max(np.abs(_tensor(wf, ctmode)-ref)) < 1e-10


def test_native_spinful_in_an_nf_and_sz_sector():
    """The ground state of this chain at Nf=4 has Sz=0, so the Nf=4, Sz=0
    sector holds the same state, and the fold there gives the same
    tensor; ctmode="full" raises, and "explicit" refuses the tuples that
    flip a spin."""
    nf = _native(dict(Nf=4))
    both = _native(dict(Nf=4, Sz=0))
    assert both.gs_energy() == pytest.approx(nf.gs_energy(), abs=1e-8)
    ref = _tensor(nf.get_gs(), "fold")
    wf = both.get_gs()
    assert np.max(np.abs(_tensor(wf)-ref)) < 1e-6
    with pytest.raises(ValueError, match="fold"):
        _tensor(wf, "full")
