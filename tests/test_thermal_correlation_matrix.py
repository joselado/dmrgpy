"""Regression coverage for Many_Body_Chain.get_correlation_matrix(mode="ED",
T=...) (entropytk/correlationentropy.py::get_correlation_matrix_finiteT).

Before the fix, this had three independent bugs:

1. Backward Boltzmann sign: `np.exp(es[i]/T)` instead of
   `np.exp(-es[i]/T)`, so the "thermal" average was actually weighted
   towards the *highest*-energy eigenstates instead of the lowest --
   confirmed by comparing against the already-correct thermal `vev`
   (post its own truncation fix, see test_thermal_vev.py), which for a
   diagonal entry should agree exactly (`cm[i,i] == vev(N[i], T=...)`).
2. A Hilbert-space-size guard based on `2**len(self.C)` (site count)
   requested unconditionally, rather than the actual Hamiltonian
   dimension against algebra.maxsize like vevtk.thermalvev.thermal_vev_ex
   does -- for n=11 sites (2**11=2048 > algebra.maxsize=2000) this
   unconditionally asked algebra.lowest_states for the *entire* spectrum
   via a sparse ARPACK call requesting k=2048 eigenvectors of a
   2048-dimensional operator, which scipy rejects outright
   ("Cannot use scipy.linalg.eig for sparse A with k >= N - 1").
3. No truncation-safety check at all (the T=0 vev bug's version at
   least printed a warning): silently wrong for any temperature/size
   combination where the default excited-state count didn't capture
   the actual Boltzmann tail.

Fixed by modeling get_correlation_matrix_finiteT directly on
thermal_vev_ex's already-correct pattern: get the true Hilbert
dimension from get_ED_obj().get_hamiltonian().shape[0], default to the
full spectrum only when that's small enough to diagonalize exactly,
correct the Boltzmann sign, and raise RuntimeError instead of silently
returning a wrong average when a truncated sum still carries
non-negligible weight in its highest-included state.
"""
import numpy as np
import pytest

from dmrgpy import fermionchain


def _chain(n, bias=0.5):
    fc = fermionchain.Fermionic_Chain(n)
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i]
    h = h + bias * fc.N[0]
    fc.set_hamiltonian(h)
    return fc


def test_correlation_matrix_diagonal_matches_thermal_vev():
    """cm[i,i] = <N_i> at finite T, which must agree with the
    independently-implemented, already-fixed thermal vev to machine
    precision -- this is exactly the sign/normalization cross-check
    that would have caught the backward-Boltzmann-sign bug (previously
    cm[0,0] and vev(N[0]) summed to ~1, i.e. exactly swapped)."""
    n = 6
    fc = _chain(n)

    for T in (0.01, 0.5, 2.0, 100.0):
        cm = fc.get_correlation_matrix(T=T, mode="ED")
        for i in range(n):
            v = fc.vev(fc.N[i], mode="ED", T=T)
            assert cm[i, i].real == pytest.approx(v.real, abs=1e-8)
            assert cm[i, i].imag == pytest.approx(0., abs=1e-8)


def test_correlation_matrix_is_hermitian():
    fc = _chain(6)
    cm = fc.get_correlation_matrix(T=1.0, mode="ED")
    assert np.allclose(cm, cm.conj().T, atol=1e-10)


def test_correlation_matrix_matches_full_spectrum_reference():
    """Independent reference built directly with numpy.linalg.eigh on
    the dense Hamiltonian (bypassing get_correlation_matrix_finiteT and
    get_correlation_matrix_zeroT entirely), not a regression-pinned
    golden value."""
    n = 6
    fc = _chain(n)
    T = 0.7

    edobj = fc.get_ED_obj()
    Hd = np.array(edobj.get_hamiltonian().todense())
    es, vs = np.linalg.eigh(Hd)
    es = es - es.min()
    P = np.exp(-es / T)
    P /= P.sum()
    Nd = [np.array(edobj.MO2matrix(fc.N[i]).todense()) for i in range(n)]

    expected_diag = np.zeros(n)
    for i in range(n):
        vals = np.array([np.conjugate(vs[:, k]) @ Nd[i] @ vs[:, k] for k in range(len(es))])
        expected_diag[i] = np.sum(P * vals).real

    cm = fc.get_correlation_matrix(T=T, mode="ED")
    assert np.diag(cm).real == pytest.approx(expected_diag, abs=1e-8)


def test_correlation_matrix_previously_crashing_large_chain_now_works():
    """n=11 sites -> Hilbert dim 2048 > algebra.maxsize=2000: the old
    `2**len(self.C)` guard unconditionally requested the full spectrum
    from a sparse ARPACK call, which scipy rejects for k>=N-1, crashing
    outright regardless of T. At low enough T the (now dimension-aware)
    default truncation is adequate and must agree with thermal vev."""
    n = 11
    fc = _chain(n)
    T = 0.01

    cm = fc.get_correlation_matrix(T=T, mode="ED")
    v = fc.vev(fc.N[0], mode="ED", T=T)
    assert cm[0, 0].real == pytest.approx(v.real, abs=1e-6)


def test_correlation_matrix_raises_instead_of_silently_truncating():
    """Explicitly requesting too few excited states for the temperature
    must raise, not silently return a wrong matrix."""
    fc = _chain(6)
    with pytest.raises(RuntimeError):
        fc.get_correlation_matrix(T=2.0, mode="ED", n=5)
