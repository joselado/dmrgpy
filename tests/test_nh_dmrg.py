"""Tests for the non-Hermitian DMRG (NH-DMRG) solver of the ITensor v3
backend (mpscpp3/chain_session.h's Chain::nhdmrg, driven by nhdmrg.py) --
a port of ITensorNHDMRG.jl's default "onesided" + "fidelity"
configuration. Each test cross-checks against exact diagonalization
(mode="ED", which diagonalizes the full non-Hermitian matrix), and one
also against the pre-existing MPS Arnoldi route
(mpsalgebra.lowest_energy_non_hermitian_arnoldi), which remains the
non-Hermitian fallback for the other backends.

Conventions checked: the targeted eigenvalue is the one with smallest
real part; (energy, psil, psir) is a biorthogonal left/right eigenpair
with <psil|psir> = 1.

NH-DMRG starts from an unseeded random MPS, so equality assertions
against ED go through moderately loose tolerances (the converged runs
observed while developing this sit at ~1e-14; 1e-6 leaves headroom
without letting a stalled run through)."""

import numpy as np
import pytest

from dmrgpy import fermionchain, spinchain, cppext

pytestmark = pytest.mark.skipif(
    not cppext.available(3),
    reason="requires the compiled mpscpp3 (ITensor v3) extension",
)


def nh_fermion_chain(n):
    """Interacting fermionic chain with a staggered imaginary potential
    (the model of examples/non_hermitian/non_hermitian_chain): complex
    spectrum, unique smallest-real-part value up to a conjugate pair."""
    fc = fermionchain.Fermionic_Chain(n)
    fc.itensor_version = 3
    fc.setup_cpp(version=3)
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i]
    for i in range(n):
        h = h + 1j * (-1)**i * fc.Cdag[i] * fc.C[i]
    for i in range(n - 1):
        h = h + (fc.N[i] - 0.5) * (fc.N[i + 1] - 0.5)
    fc.set_hamiltonian(h)
    return fc, h


def nh_pt_spin_chain(n, g=0.3):
    """PT-symmetric XX chain with a staggered imaginary field (the model
    class of Yamamoto et al., PRB 105, 205125): several eigenvalues share
    the smallest real part (a complex-conjugate pair plus real ones), the
    degenerate case that requires the left solve to be anchored to the
    right solve's eigenvalue (see arnoldi_smallest_real's Sel comment in
    mpscpp3/chain_session.h)."""
    sc = spinchain.Spin_Chain(["S=1/2"] * n)
    sc.itensor_version = 3
    sc.setup_cpp(version=3)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1]
    for i in range(n):
        h = h + 1j * g * (-1)**i * sc.Sz[i]
    sc.set_hamiltonian(h)
    return sc, h


def test_nhdmrg_matches_ed_fermionic_chain():
    fc, h = nh_fermion_chain(4)
    es_ed = fc.get_excited(mode="ED", n=4)
    e, psil, psir = fc.nhdmrg()
    # smallest real part reproduced (the ED list is sorted by real part)
    assert e.real == pytest.approx(es_ed[0].real, abs=1e-6)
    # and the energy is an actual eigenvalue (either member of the
    # conjugate pair that shares the smallest real part is acceptable)
    assert min(abs(e - x) for x in es_ed) == pytest.approx(0.0, abs=1e-6)


def test_nhdmrg_left_right_eigenpair():
    fc, h = nh_fermion_chain(4)
    e, psil, psir = fc.nhdmrg()
    # biorthonormal pair
    assert psil.dot(psir) == pytest.approx(1.0, abs=1e-8)
    # right eigenvector: H|r> = E|r>
    r = h * psir - e * psir
    assert abs(r.dot(r))**0.5 == pytest.approx(0.0, abs=1e-6)
    # left eigenvector: Hdag|l> = conj(E)|l>
    l = h.get_dagger() * psil - np.conj(e) * psil
    assert abs(l.dot(l))**0.5 == pytest.approx(0.0, abs=1e-6)
    # biorthogonal Rayleigh quotient reproduces the eigenvalue
    assert psil.aMb(h, psir) / psil.dot(psir) == pytest.approx(e, abs=1e-8)


def test_nhdmrg_matches_arnoldi():
    """NH-DMRG and the pre-existing MPS Arnoldi route agree on the same
    smallest-real-part eigenvalue (Arnoldi converges far less tightly
    than NH-DMRG -- ~1e-3 vs ~1e-14 observed on this model -- so the
    cross-tolerance is Arnoldi's, and each is also pinned to ED at its
    own accuracy)."""
    from dmrgpy import mpsalgebra
    fc, h = nh_fermion_chain(4)
    es_ed = fc.get_excited(mode="ED", n=4)
    e, psil, psir = fc.nhdmrg()
    ea, wfa = mpsalgebra.lowest_energy_non_hermitian_arnoldi(fc, h, n=1)
    assert e.real == pytest.approx(ea[0].real, abs=5e-2)
    assert min(abs(ea[0] - x) for x in es_ed) == pytest.approx(0.0, abs=5e-2)
    assert min(abs(e - x) for x in es_ed) == pytest.approx(0.0, abs=1e-6)


def test_nhdmrg_pt_symmetric_spin_chain():
    sc, h = nh_pt_spin_chain(6)
    es_ed = sc.get_excited(mode="ED", n=6)
    e, psil, psir = sc.nhdmrg()
    assert e.real == pytest.approx(es_ed[0].real, abs=1e-6)
    assert min(abs(e - x) for x in es_ed) == pytest.approx(0.0, abs=1e-6)
    # converged to a true eigenpair, not a variational stall (the
    # non-Hermitian "energy" is not a variational bound, so the residual
    # is the meaningful convergence certificate)
    r = h * psir - e * psir
    assert abs(r.dot(r))**0.5 == pytest.approx(0.0, abs=1e-6)


def test_gs_energy_routes_to_nhdmrg_on_v3():
    """For a non-Hermitian Hamiltonian on itensor_version=3, gs_energy()
    now runs NH-DMRG (groundstate.py's non-Hermitian branch) instead of
    the Arnoldi route, stores the right eigenvector as wf0, and returns
    the smallest-real-part eigenvalue."""
    fc, h = nh_fermion_chain(4)
    es_ed = fc.get_excited(mode="ED", n=4)
    e0 = fc.gs_energy()
    assert e0.real == pytest.approx(es_ed[0].real, abs=1e-6)
    assert fc.computed_gs
    # wf0 holds the right eigenvector
    r = h * fc.wf0 - e0 * fc.wf0
    assert abs(r.dot(r))**0.5 == pytest.approx(0.0, abs=1e-6)
