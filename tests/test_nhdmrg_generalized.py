"""Regression coverage for the non-Hermitian generalized-eigenvalue
NH-DMRG solver: H|psi_R>=lambda*A|psi_R> for a possibly non-Hermitian H
and a Hermitian positive-definite metric operator A -- the non-Hermitian
counterpart of tests/test_dmrg_generalized.py, generalizing NH-DMRG
(tests/test_nh_dmrg.py) exactly the way gs_energy_generalized()
generalizes plain gs_energy(). See pyitensor/nhdmrg.py's
nhdmrg_generalized() and mpscpp3/chain_session.h's
Chain::nhdmrg_generalized for the self-consistent Lagrange-multiplier
algorithm (complex lambda, biorthogonal Rayleigh quotient) -- the same
algorithm on both backends, so every accuracy test below is parametrized
over itensor_version in ("python", 3); the v3 case is skipped
automatically when the compiled extension isn't available. mpscpp2
(itensor_version=2) has no analogous session method yet.

Ground truth is scipy.linalg.eig's general (non-symmetric) generalized
eigensolver, built from the ED sparse matrices of H and A -- unlike
test_dmrg_generalized.py's scipy.linalg.eigh (Hermitian-definite only),
since H is non-Hermitian here. NH-DMRG starts from an unseeded random
MPS (same as plain nhdmrg(), see test_nh_dmrg.py's own note), so
equality assertions against ED use moderately loose tolerances.
"""
import numpy as np
import pytest
import scipy.linalg as sla

from dmrgpy import cppext, fermionchain
from dmrgpy.multioperator import identity as mo_identity
from dmrgpy.nhdmrg import nhdmrg_generalized as _nhdmrg_generalized_toplevel

ITENSOR_VERSIONS = [
    "python",
    pytest.param(3, marks=pytest.mark.skipif(
        not cppext.available(3), reason="ITensor v3 extension not compiled")),
]


def _nh_generalized_fermion_problem(n=4, seed=2):
    """H a non-Hermitian interacting fermion chain (staggered imaginary
    on-site potential, the model of examples/non_hermitian/
    non_hermitian_chain, same as test_nh_dmrg.py's own nh_fermion_chain),
    A = 1 + 0.5*sum(N_i) a Hermitian positive-definite metric."""
    rng = np.random.RandomState(seed)
    fc = fermionchain.Fermionic_Chain(n)
    h = 0
    for i in range(n - 1):
        h = h + rng.random() * (fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i])
    for i in range(n):
        h = h + 1j * (-1) ** i * 0.3 * fc.N[i]
    for i in range(n - 1):
        h = h + rng.random() * (fc.N[i] - 0.5) * (fc.N[i + 1] - 0.5)
    m = 1
    for i in range(n):
        m = m + 0.5 * fc.N[i]
    fc.set_hamiltonian(h)
    return fc, h, m


def _nh_generalized_ground_truth(fc, h, m):
    Hmat = fc.get_ED_obj().MO2matrix(h)
    Mmat = fc.get_ED_obj().MO2matrix(m)
    Hmat = Hmat.toarray() if hasattr(Hmat, "toarray") else np.array(Hmat)
    Mmat = Mmat.toarray() if hasattr(Mmat, "toarray") else np.array(Mmat)
    w = sla.eig(Hmat, Mmat, right=False)
    return w[np.argsort(w.real)]


def _setup(fc, itensor_version, maxm=40, nsweeps=14, cutoff=1e-12):
    if itensor_version == "python":
        fc.setup_python()
    else:
        fc.setup_cpp(version=itensor_version)
    fc.maxm = maxm
    fc.nsweeps = nsweeps
    fc.cutoff = cutoff
    return fc


@pytest.mark.parametrize("itensor_version", ITENSOR_VERSIONS)
def test_gs_energy_generalized_nonhermitian_matches_ed(itensor_version):
    """The non-Hermitian generalized NH-DMRG eigenvalue (smallest real
    part) must agree with scipy.linalg.eig's general generalized
    eigensolver."""
    fc, h, m = _nh_generalized_fermion_problem(n=4, seed=2)
    w = _nh_generalized_ground_truth(fc, h, m)
    _setup(fc, itensor_version)
    lam = fc.gs_energy_generalized(m)
    assert abs(lam - w[0]) < 1e-3


@pytest.mark.parametrize("itensor_version", ITENSOR_VERSIONS)
def test_nhdmrg_generalized_matches_ed_and_is_biorthogonal(itensor_version):
    """Drive the top-level nhdmrg_generalized() directly (not through
    gs_energy_generalized()'s retry wrapper) to check both the eigenvalue
    and the biorthogonal normalization <psi_L|psi_R>=1 of the returned
    pair, using the same backend-agnostic mps.MPS.dot() every backend
    supports (mirrors test_nh_dmrg.py's own convention)."""
    fc, h, m = _nh_generalized_fermion_problem(n=4, seed=2)
    w = _nh_generalized_ground_truth(fc, h, m)
    _setup(fc, itensor_version)
    lam, psil, psir = _nhdmrg_generalized_toplevel(fc, m, H=h,
            krylovdim=20, restarts=3, ntries=3)
    assert abs(lam - w[0]) < 1e-4
    assert abs(psil.dot(psir) - 1.0) < 1e-6


@pytest.mark.parametrize("itensor_version", ITENSOR_VERSIONS)
def test_nhdmrg_generalized_residual_is_small(itensor_version):
    """A genuine eigen-residual certificate, not just agreement with ED:
    ||H|psi_R>-lambda*A|psi_R>|| must be small for the returned pair."""
    fc, h, m = _nh_generalized_fermion_problem(n=4, seed=2)
    _setup(fc, itensor_version)
    lam, psil, psir = _nhdmrg_generalized_toplevel(fc, m, H=h,
            krylovdim=20, restarts=3, ntries=3)
    r = h * psir - lam * (m * psir)
    resid_norm = abs(r.dot(r)) ** 0.5
    assert resid_norm < 1e-3


@pytest.mark.parametrize("itensor_version", ITENSOR_VERSIONS)
def test_nhdmrg_generalized_identity_matches_plain_nhdmrg(itensor_version):
    """A = identity must reduce the non-Hermitian generalized problem
    exactly to plain nhdmrg()'s own eigenvalue (nhdmrg_generalized()'s
    own docstring: A=Id makes lambda the ordinary biorthogonal Rayleigh
    quotient, i.e. plain NH-DMRG).

    This test problem's spectrum has a complex-conjugate-degenerate pair
    sharing the smallest real part (confirmed directly: repeated trials
    reproducibly landed on either -1.208732-0.348597j or its conjugate,
    with real parts always agreeing to ~1e-15) -- the same structural
    ambiguity test_nh_dmrg.py's own nh_fermion_chain tests explicitly
    account for ("either member of the conjugate pair ... is
    acceptable"), since each independently-random-seeded NH-DMRG run
    (both nhdmrg() and nhdmrg_generalized() start from an unseeded random
    MPS) may land on either member with no way to control which. Compare
    against whichever of e_plain/conj(e_plain) is closer, not e_plain
    alone."""
    fc_plain, h_plain, _ = _nh_generalized_fermion_problem(n=4, seed=3)
    _setup(fc_plain, itensor_version)
    e_plain, _, _ = fc_plain.nhdmrg()

    fc_gen, h_gen, _ = _nh_generalized_fermion_problem(n=4, seed=3)
    _setup(fc_gen, itensor_version)
    lam = fc_gen.gs_energy_generalized(mo_identity())
    err = min(abs(lam - e_plain), abs(lam - np.conj(e_plain)))
    assert err < 1e-3


@pytest.mark.parametrize("itensor_version", ITENSOR_VERSIONS)
def test_gs_energy_generalized_nonhermitian_leaves_computed_gs_true(itensor_version):
    """Same regression as test_dmrg_generalized.py's Hermitian-path
    version: computed_gs must end up True (not silently False) after the
    call, and a subsequent gs_energy() must return the cached value
    rather than quietly re-solving to something else."""
    fc, h, m = _nh_generalized_fermion_problem(n=4, seed=2)
    _setup(fc, itensor_version)
    lam = fc.gs_energy_generalized(m)
    assert fc.computed_gs is True
    assert fc.gs_energy() == lam


def test_gs_energy_generalized_nonhermitian_requires_supported_backend():
    """The non-Hermitian generalized solver exists on itensor_version
    "python" and 3 (like the Hermitian path) but not mpscpp2
    (itensor_version=2), which has no Chain.nhdmrg_generalized yet."""
    fc, h, m = _nh_generalized_fermion_problem(n=4, seed=2)
    fc.setup_cpp(version=2)
    with pytest.raises(NotImplementedError):
        fc.gs_energy_generalized(m)


def test_nhdmrg_generalized_v3_short_chain_does_not_crash():
    """Unlike the Hermitian path (which calls ITensor v3's own dmrg() and
    aborts the whole process for chains shorter than 3 sites, see
    mode.py's own guard), NH-DMRG's two-site sweep is hand-rolled
    directly against arnoldi_smallest_real/manual ITensor contractions
    and never calls dmrg() at all -- so a 2-site chain must run this
    method to completion instead of being rejected or crashing."""
    if not cppext.available(3):
        pytest.skip("ITensor v3 extension not compiled")
    fc, h, m = _nh_generalized_fermion_problem(n=2, seed=2)
    _setup(fc, 3, maxm=20, nsweeps=8)
    lam = fc.gs_energy_generalized(m)  # must not raise or crash
    assert np.isfinite(lam.real) and np.isfinite(lam.imag)


def test_nhdmrg_generalized_lam0_bare_float_nan_is_treated_as_unset():
    """Regression test for a code-review finding: nhdmrg_generalized()'s
    original lam0 "unset" check was `isinstance(lam, complex) and
    np.isnan(lam.real)`, so a caller passing a *bare float* NaN (a
    reasonable thing to do given this function's own docstring wording,
    "None or a NaN real part") failed the isinstance(complex) gate and
    fed a live NaN straight into the first outer sweep instead of
    triggering the documented auto-seed fallback. lam0_is_unset()
    (dmrg.py, shared with the Hermitian solver) now coerces via
    complex(lam0) instead, so a bare float NaN must be treated exactly
    like None."""
    fc, h, m = _nh_generalized_fermion_problem(n=4, seed=2)
    w = _nh_generalized_ground_truth(fc, h, m)
    _setup(fc, "python")
    lam = fc.gs_energy_generalized(m, lam0=float('nan'))
    assert abs(lam - w[0]) < 1e-3


def test_nhdmrg_generalized_retry_loop_survives_a_transient_runtime_error(monkeypatch):
    """Regression test for a code-review finding: the ntries retry loop
    had no try/except around the session call, so the new "<psi_L|A|
    psi_R> collapsed to ~0" guard (raised on an unlucky random draw)
    aborted every remaining attempt instead of being retried like an
    ordinary resid>=tol failure -- defeating the whole point of
    ntries>1 (each attempt starts from an independent fresh random
    state). Simulated here via monkeypatch (reliably reproducing a real
    near-null-space draw is not deterministic) rather than trying to
    engineer a genuinely bad random seed."""
    fc, h, m = _nh_generalized_fermion_problem(n=4, seed=2)
    _setup(fc, "python")
    real_method = fc._session.nhdmrg_generalized
    calls = {"n": 0}

    def flaky(*args, **kwargs):
        calls["n"] += 1
        if calls["n"] == 1:
            raise RuntimeError("simulated near-null-space collapse")
        return real_method(*args, **kwargs)

    monkeypatch.setattr(fc._session, "nhdmrg_generalized", flaky)
    lam, psil, psir = _nhdmrg_generalized_toplevel(fc, m, H=h, ntries=3)
    assert calls["n"] >= 2  # the first (failing) attempt was retried, not fatal
    assert abs(psil.dot(psir) - 1.0) < 1e-6


def test_nhdmrg_generalized_all_attempts_failing_raises_clear_error(monkeypatch):
    """If every attempt (up to ntries) hits the collapse guard, the
    top-level wrapper must raise a clear RuntimeError of its own rather
    than crashing on an unbound `best` variable or silently returning
    None-derived garbage."""
    fc, h, m = _nh_generalized_fermion_problem(n=4, seed=2)
    _setup(fc, "python")

    def always_fails(*args, **kwargs):
        raise RuntimeError("simulated near-null-space collapse")

    monkeypatch.setattr(fc._session, "nhdmrg_generalized", always_fails)
    with pytest.raises(RuntimeError):
        _nhdmrg_generalized_toplevel(fc, m, H=h, ntries=3)
