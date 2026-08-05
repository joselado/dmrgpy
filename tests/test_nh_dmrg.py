"""Tests for the non-Hermitian DMRG (NH-DMRG) solver -- a port of
ITensorNHDMRG.jl's default "onesided" + "fidelity" configuration,
implemented on all three session backends: mpscpp3 (the annotated
original, mpscpp3/chain_session.h's Chain::nhdmrg), mpscpp2 (its v2-API
back-port), and pyitensor (pyitensor/nhdmrg.py). itensor_version=
"julia_live" is covered too, but is not a port at all -- it calls the
real ITensorNHDMRG.jl package, with two reconciling fixes on top (a
left-vector conjugation and a real-part-degeneracy tie-break, see
mpsjulialive/nhdmrg.jl's header and nh_asymmetric_hopping_chain's
docstring below). All four are driven by the shared
retry/certificate wrapper in nhdmrg.py. Each test cross-checks against
exact diagonalization (mode="ED", which diagonalizes the full
non-Hermitian matrix), and one also against the pre-existing MPS Arnoldi
route (mpsalgebra.lowest_energy_non_hermitian_arnoldi), the previous
non-Hermitian fallback.

Conventions checked: the targeted eigenvalue is the one with smallest
real part; (energy, psil, psir) is a biorthogonal left/right eigenpair
with <psil|psir> = 1.

NH-DMRG starts from an unseeded random MPS, so equality assertions
against ED go through moderately loose tolerances; residual assertions
use 1e-3, comfortably above the driver's own acceptance certificate
(relative 1e-4 -- an accepted run can sit just below that) and far below
a stalled run's ~1e-1."""

import numpy as np
import pytest

from dmrgpy import fermionchain, spinchain, cppext

from _helpers import julia_live_param, setup_backend

VERSIONS = [
    pytest.param(2, marks=pytest.mark.skipif(
        not cppext.available(2),
        reason="requires the compiled mpscpp2 (ITensor v2) extension")),
    pytest.param(3, marks=pytest.mark.skipif(
        not cppext.available(3),
        reason="requires the compiled mpscpp3 (ITensor v3) extension")),
    pytest.param("python", id="python"),
    julia_live_param(),
]


def nh_fermion_chain(n, version):
    """Interacting fermionic chain with a staggered imaginary potential
    (the model of examples/non_hermitian/non_hermitian_chain): complex
    spectrum, unique smallest-real-part value up to a conjugate pair."""
    fc = fermionchain.Fermionic_Chain(n)
    setup_backend(fc, version)
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i]
    for i in range(n):
        h = h + 1j * (-1)**i * fc.Cdag[i] * fc.C[i]
    for i in range(n - 1):
        h = h + (fc.N[i] - 0.5) * (fc.N[i + 1] - 0.5)
    fc.set_hamiltonian(h)
    return fc, h


def nh_pt_spin_chain(n, version, g=0.3):
    """PT-symmetric XX chain with a staggered imaginary field (the model
    class of Yamamoto et al., PRB 105, 205125): several eigenvalues share
    the smallest real part (a complex-conjugate pair plus real ones), the
    degenerate case that requires the left solve to be anchored to the
    right solve's eigenvalue (see arnoldi_smallest_real's Sel comment in
    mpscpp3/chain_session.h)."""
    sc = spinchain.Spin_Chain(["S=1/2"] * n)
    setup_backend(sc, version)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1]
    for i in range(n):
        h = h + 1j * g * (-1)**i * sc.Sz[i]
    sc.set_hamiltonian(h)
    return sc, h


def nh_asymmetric_hopping_chain(n, version, tr=1.0, tl=0.6, gamma=0.7):
    """Hatano-Nelson-style chain: *asymmetric* hopping (tr != tl) plus a
    staggered imaginary on-site potential.

    The point of this model is that it is genuinely non-symmetric --
    H^T != H, on top of H^dagger != H -- while still having a
    complex-conjugate pair tied for the smallest real part. Both models
    above are complex *symmetric* (H^T == H: their hoppings are symmetric
    and every non-Hermitian piece is diagonal), and for such a
    Hamiltonian the transpose left eigenvector and the adjoint one
    coincide up to a complex conjugation, so a solver that returns the
    wrong one of the two still passes every assertion made about them.
    That is not hypothetical: itensor_version="julia_live" calls the real
    ITensorNHDMRG.jl package, and both of mpsjulialive/nhdmrg.jl's
    julia_live-specific fixes are invisible without a model like this one
    (see that file's header for both in full):

    - its "adjoint" sweep is against swapprime(H,0=>1) == H^T rather than
      H^dagger, so the returned left vector needs a conjugation
      (nh_biorthogonal_pair); without it the left residual sits at ~2
      while the right one is ~1e-15;
    - nothing in it ties the left solve to whichever member of the
      real-part-degenerate pair the right solve picked, so the two
      converged to *different* eigenvalues -- deterministically, every
      attempt from every random start, with overlap exactly 0
      (nhdmrg_solve's exp(i*theta) tie-break).

    tl=0.6 specifically: confirmed to be one of the values where the
    second failure mode fires on every attempt without the tie-break, so
    this is a real regression guard rather than a probabilistic one.
    """
    fc = fermionchain.Fermionic_Chain(n)
    setup_backend(fc, version)
    h = 0
    for i in range(n - 1):
        h = h + tr * fc.Cdag[i] * fc.C[i + 1] + tl * fc.Cdag[i + 1] * fc.C[i]
    for i in range(n):
        h = h + 1j * gamma * (-1)**i * fc.Cdag[i] * fc.C[i]
    for i in range(n - 1):
        h = h + 0.5 * (fc.N[i] - 0.5) * (fc.N[i + 1] - 0.5)
    fc.set_hamiltonian(h)
    return fc, h


@pytest.mark.parametrize("version", VERSIONS)
def test_nhdmrg_asymmetric_hopping_left_right_eigenpair(version):
    """The left eigenvector must satisfy the *adjoint* equation
    H^dagger|psil> = conj(E)|psil>, not the transpose one -- checked on a
    Hamiltonian where the two genuinely differ (see
    nh_asymmetric_hopping_chain's docstring)."""
    fc, h = nh_asymmetric_hopping_chain(4, version)
    es_ed = fc.get_excited(mode="ED", n=4)
    e, psil, psir = fc.nhdmrg()
    assert e.real == pytest.approx(es_ed[0].real, abs=1e-6)
    assert min(abs(e - x) for x in es_ed) == pytest.approx(0.0, abs=1e-6)
    assert psil.dot(psir) == pytest.approx(1.0, abs=1e-8)
    r = h * psir - e * psir
    assert abs(r.dot(r))**0.5 == pytest.approx(0.0, abs=1e-3)
    l = h.get_dagger() * psil - np.conj(e) * psil
    assert abs(l.dot(l))**0.5 == pytest.approx(0.0, abs=1e-3)


def test_nhdmrg_requires_a_usable_session():
    """itensor_version 2/3 with no compiled extension (self._session left
    as None by sites.py's initialize()) must raise a clear RuntimeError,
    not an AttributeError from calling into a None session -- the same
    guard gs_energy_generalized() and nhdmrg_generalized() already had,
    which nhdmrg() was missing."""
    fc = fermionchain.Fermionic_Chain(4)
    h = fc.Cdag[0] * fc.C[1] + fc.Cdag[1] * fc.C[0] + 1j * fc.N[0]
    fc.set_hamiltonian(h)
    fc.setup_cpp(version=3)
    fc._session = None  # simulate "extension not compiled" regardless of env
    with pytest.raises(RuntimeError):
        fc.nhdmrg()


def test_nhdmrg_all_attempts_failing_reports_the_underlying_error(monkeypatch):
    """nhdmrg()'s retry loop swallows RuntimeError to redraw a bad random
    start, so a *deterministic* failure gets retried ntries times and then
    reported by the loop itself. That report must carry the attempt's own
    message (and __cause__) through: every backend raises RuntimeError for
    its own reasons -- an ITensor Error() surfaces as one through pybind11
    -- and stating only this driver's guess would send the user tuning
    nsweeps/maxm/krylovdim for a problem none of them affect."""
    fc, h = nh_fermion_chain(4, "python")

    def always_fails(*args, **kwargs):
        raise RuntimeError("simulated deterministic backend failure")

    monkeypatch.setattr(fc._session, "nhdmrg", always_fails)
    with pytest.raises(RuntimeError) as excinfo:
        fc.nhdmrg(ntries=2)
    assert "simulated deterministic backend failure" in str(excinfo.value)
    assert excinfo.value.__cause__ is not None


@pytest.mark.parametrize("version", VERSIONS)
def test_nhdmrg_matches_ed_fermionic_chain(version):
    fc, h = nh_fermion_chain(4, version)
    es_ed = fc.get_excited(mode="ED", n=4)
    e, psil, psir = fc.nhdmrg()
    # smallest real part reproduced (the ED list is sorted by real part)
    assert e.real == pytest.approx(es_ed[0].real, abs=1e-6)
    # and the energy is an actual eigenvalue (either member of the
    # conjugate pair that shares the smallest real part is acceptable)
    assert min(abs(e - x) for x in es_ed) == pytest.approx(0.0, abs=1e-6)


@pytest.mark.parametrize("version", VERSIONS)
def test_nhdmrg_left_right_eigenpair(version):
    fc, h = nh_fermion_chain(4, version)
    e, psil, psir = fc.nhdmrg()
    # biorthonormal pair
    assert psil.dot(psir) == pytest.approx(1.0, abs=1e-8)
    # right eigenvector: H|r> = E|r>
    r = h * psir - e * psir
    assert abs(r.dot(r))**0.5 == pytest.approx(0.0, abs=1e-3)
    # left eigenvector: Hdag|l> = conj(E)|l>
    l = h.get_dagger() * psil - np.conj(e) * psil
    assert abs(l.dot(l))**0.5 == pytest.approx(0.0, abs=1e-3)
    # biorthogonal Rayleigh quotient reproduces the eigenvalue
    assert psil.aMb(h, psir) / psil.dot(psir) == pytest.approx(e, abs=1e-8)


def test_nhdmrg_matches_arnoldi():
    """NH-DMRG and the pre-existing MPS Arnoldi route agree on the same
    smallest-real-part eigenvalue (Arnoldi converges far less tightly
    than NH-DMRG -- ~1e-3 vs ~1e-14 observed on this model -- so the
    cross-tolerance is Arnoldi's, and each is also pinned to ED at its
    own accuracy)."""
    from dmrgpy import mpsalgebra
    if not cppext.available(3):
        pytest.skip("requires the compiled mpscpp3 (ITensor v3) extension")
    fc, h = nh_fermion_chain(4, 3)
    es_ed = fc.get_excited(mode="ED", n=4)
    e, psil, psir = fc.nhdmrg()
    ea, wfa = mpsalgebra.lowest_energy_non_hermitian_arnoldi(fc, h, n=1)
    assert e.real == pytest.approx(ea[0].real, abs=5e-2)
    assert min(abs(ea[0] - x) for x in es_ed) == pytest.approx(0.0, abs=5e-2)
    assert min(abs(e - x) for x in es_ed) == pytest.approx(0.0, abs=1e-6)


@pytest.mark.parametrize("version", VERSIONS)
def test_nhdmrg_pt_symmetric_spin_chain(version):
    sc, h = nh_pt_spin_chain(6, version)
    es_ed = sc.get_excited(mode="ED", n=6)
    e, psil, psir = sc.nhdmrg()
    assert e.real == pytest.approx(es_ed[0].real, abs=1e-6)
    assert min(abs(e - x) for x in es_ed) == pytest.approx(0.0, abs=1e-6)
    # converged to a true eigenpair, not a variational stall (the
    # non-Hermitian "energy" is not a variational bound, so the residual
    # is the meaningful convergence certificate)
    r = h * psir - e * psir
    assert abs(r.dot(r))**0.5 == pytest.approx(0.0, abs=1e-3)


@pytest.mark.parametrize("version", VERSIONS)
def test_gs_energy_routes_to_nhdmrg(version):
    """For a non-Hermitian Hamiltonian, gs_energy() runs NH-DMRG
    (groundstate.py's non-Hermitian branch) instead of the Arnoldi route,
    stores the right eigenvector as wf0, and returns the
    smallest-real-part eigenvalue. On the session backends (v2, v3, pure
    Python) that means nhdmrg.py's driver; julia_live reaches the same
    ITensorNHDMRG.jl solver through its own, older branch in
    groundstate.py (mpsjulialive/groundstate.py's get_gs_dmrg with
    ishermitian=False), which keeps only the right eigenvector -- the
    biorthogonal pair is what chain.nhdmrg() is for there."""
    fc, h = nh_fermion_chain(4, version)
    es_ed = fc.get_excited(mode="ED", n=4)
    e0 = fc.gs_energy()
    assert e0.real == pytest.approx(es_ed[0].real, abs=1e-6)
    assert fc.computed_gs
    # wf0 holds the right eigenvector
    r = h * fc.wf0 - e0 * fc.wf0
    assert abs(r.dot(r))**0.5 == pytest.approx(0.0, abs=1e-3)
