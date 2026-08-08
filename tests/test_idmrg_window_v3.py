"""Coverage for the native ITensor v3 (mpscpp3) port of the IBC-window
real-time dynamical correlator (Milsted/Vanderstraeten, arXiv:1804.09163,
Sec. V.1) -- `Chain::td_dynamical_correlator_window` (mpscpp3/chain_session.h)
and `Infinite_Many_Body_Chain.td_dynamical_correlator(itensor_version=3)`
(infinitechain.py). See `tests/test_idmrg_window.py` for the original
pyitensor-only coverage this mirrors, and `pyitensor/idmrg_window.py`'s own
module docstring for the algorithm.

Model choice: pure isotropic S=1/2 Heisenberg (n_uc=1, no onsite field
term), matching `examples/idmrg/heisenberg_infinite_python_VS_v3/main.py`
(the existing, known-good python-vs-v3 `idmrg_ground_state` regression
check) -- deliberately *not* the transverse-field model
`test_idmrg_window.py` itself uses (`h_intra=[[1.4,["Sz",0]]]`): confirmed
directly (reproducible on an unmodified checkout, unrelated to this PR)
that `Chain::idmrg_ground_state` diverges (energy density growing without
bound every macro-iteration) for that same onsite-field model on the v3
backend specifically -- a genuine, pre-existing v3-only bug in its onsite
("IdmrgOnsite") term handling, not something this module's own window/TDVP
code triggers or fixes. Flagged as a known follow-up, out of scope here.

Since both backends seed their local two-site solve from an unseeded
random seed (see idmrg.py's own `_local_two_site_solve` and
`arnoldi_smallest_real`'s own random start), two independent
`itensor_version="python"` and `itensor_version=3` runs converge to
independently-gauged realizations of the same physical ground state, not
bit-identical ones -- cross-backend checks below use loose tolerances
(matching this test suite's own established convention for iDMRG,
e.g. `test_idmrg_window.py`'s own `abs=0.15` energy-density checks) rather
than exact agreement, and are run at well-converged (maxm=30) settings so
run-to-run gauge/seed variation stays small.
"""
import numpy as np
import pytest

from dmrgpy import cppext
from dmrgpy.infinitechain import Infinite_Spin_Chain

pytestmark = pytest.mark.skipif(
    not cppext.available(3), reason="ITensor v3 extension not compiled")

MAXM, CUTOFF, MAXITER, ETOL, NITER, RESTARTS = 30, 1e-12, 60, 1e-9, 30, 2


def _build(itensor_version):
    """Well-converged (slow, ~minutes) -- only for tests that actually
    compare physics between the two backends and need run-to-run gauge/
    seed variation to stay small (see this module's own docstring)."""
    c = Infinite_Spin_Chain(["1/2"], itensor_version=itensor_version)
    c.gs_method = "idmrg"  # td_dynamical_correlator needs the growing
                           # algorithm's own converged-environment
                           # snapshot -- not built by gs_method="vumps"
                           # (the default since 2026-08-08)
    c.maxm, c.cutoff, c.maxiter = MAXM, CUTOFF, MAXITER
    c.etol, c.niter, c.restarts = ETOL, NITER, RESTARTS
    h = c.SxC[0] * c.SxR[0] + c.SyC[0] * c.SyR[0] + c.SzC[0] * c.SzR[0]
    c.set_hamiltonian(h)
    c.gs_energy()
    return c


def _build_fast(itensor_version):
    """Cheap, loosely-converged -- for tests that only need *some* valid
    converged environment snapshot to exercise plumbing (error handling,
    dispatch/shape), not physical accuracy."""
    c = Infinite_Spin_Chain(["1/2"], itensor_version=itensor_version)
    c.gs_method = "idmrg"  # see _build's own comment
    c.maxm, c.cutoff, c.maxiter = 8, 1e-10, 30
    c.etol, c.niter, c.restarts = 1e-6, 30, 1
    h = c.SxC[0] * c.SxR[0] + c.SyC[0] * c.SyR[0] + c.SzC[0] * c.SzR[0]
    c.set_hamiltonian(h)
    c.gs_energy()
    return c


def test_v3_idmrg_ground_state_matches_python_energy_density():
    """Sanity check that this test's own model/parameters converge
    consistently on both backends before trusting any td_dynamical_correlator
    comparison built on top -- mirrors
    heisenberg_infinite_python_VS_v3/main.py's own 1e-6 tolerance."""
    cp = _build("python")
    cv = _build(3)
    assert cv.e0 == pytest.approx(cp.e0, abs=1e-5)


def test_td_dynamical_correlator_v3_requires_gs_energy_first():
    """Chain::td_dynamical_correlator_window raises a clear RuntimeError
    (not a crash/AttributeError) when called before idmrg_ground_state --
    exercised directly against a fresh mpscpp3 Chain (bypassing
    infinitechain.py's own auto-gs_energy() convenience, which would
    otherwise always satisfy this precondition first)."""
    backend = cppext.get_backend(3)
    chain = backend.Chain([2])  # site_types: one spin-1/2, matches SpinX's convention
    with pytest.raises(RuntimeError):
        chain.td_dynamical_correlator_window(
            2, "Sz", "Sz", 0.1, 2, [0], 20, 1e-10, 50, True, 0)


def test_td_dynamical_correlator_v3_rejects_x_beyond_window():
    """x_values must keep center+x within the window's own explicit range
    -- this port doesn't implement idmrg_window.py's own beyond-window
    padding (see td_dynamical_correlator_window's own docstring), so an
    out-of-range x must raise a clear error rather than reading
    out-of-bounds/silently truncating."""
    cv = _build_fast(3)
    n_window = 2  # 2 sites total for n_uc=1 -- x=+5 is far outside
    with pytest.raises(RuntimeError):
        cv._session3.td_dynamical_correlator_window(
            n_window, "Sz", "Sz", 0.1, 1, [5], 20, 1e-10, 50, True, 0)


def test_td_dynamical_correlator_v3_eshift_insensitive_to_evolution_maxdim():
    """Basic sanity/consistency coverage for a real (now fixed) issue
    found by code review in the `eshift` fix itself (not the original
    pre-fix bug): the throwaway t=0 `tdvp()` call that measures `eshift`
    (`Chain::td_dynamical_correlator_window`, `mpscpp3/chain_session.h`)
    originally reused the *caller's* own `maxdim`/`cutoff` for its own
    measurement sweep -- but TDVPWorker's per-bond SVD split still runs
    (and truncates) even at `t=0` (`exp(0*Heff)=Id` makes the *local
    update* a no-op, not the split), so a caller-chosen `maxdim` smaller
    than the window's own already-converged bond dimension would silently
    truncate the throwaway copy *before* its energy is read off, biasing
    `eshift`. Fixed by measuring `eshift` with a fixed, generous
    `maxdim`/`cutoff=0` regardless of the caller's own evolution settings
    (mirroring pyitensor's own `window_total_energy`, which contracts
    exactly with no truncation at all).

    NOT a strong regression test for this specific fix: confirmed
    directly (temporarily reintroducing the bug and rebuilding) that at
    `_build_fast`'s own modest `maxm=8`, the resulting eshift bias is too
    small relative to ordinary evolution-truncation noise at `maxdim=4`
    to reliably fail a loose tolerance here -- isolating it cleanly would
    need a larger `maxm` and a more extreme `maxdim` mismatch, adding
    real runtime to an already-slow test module (see this module's own
    docstring), judged disproportionate for now. What this test does
    check: `td_dynamical_correlator_window` stays finite and reasonably
    self-consistent across different evolution `maxdim` choices on the
    same converged chain -- a real, if modest, sanity check."""
    n_window, dt, nt = 4, 0.1, 3
    cv = _build_fast(3)  # maxm=8

    ts_lo, xs_lo, S_lo = cv._session3.td_dynamical_correlator_window(
        n_window, "Sz", "Sz", dt, nt, [0, 1], 4, 1e-10, 50, True, 0)
    ts_hi, xs_hi, S_hi = cv._session3.td_dynamical_correlator_window(
        n_window, "Sz", "Sz", dt, nt, [0, 1], 30, 1e-10, 50, True, 0)
    S_lo, S_hi = np.array(S_lo), np.array(S_hi)

    assert np.all(np.isfinite(S_lo)) and np.all(np.isfinite(S_hi))
    assert np.max(np.abs(S_lo - S_hi)) < 0.1


def test_td_dynamical_correlator_v3_reproducible_across_convergence():
    """`Chain::td_dynamical_correlator_window`'s own `idmrg_HL_`/`idmrg_HR_`
    are, like their pyitensor counterparts, not energy-baseline-subtracted
    -- pre-fix, TDVP-evolving the (perturbed) window under them directly
    carried a spurious, run-dependent global phase (see
    tests/test_idmrg_window_free_fermion.py's own module docstring for the
    full derivation/pyitensor-side confirmation of this same bug). Fixed
    here with a post-hoc `exp(+i*eshift*t)` correction (`eshift` measured
    via a throwaway t=0 `tdvp()` call, see td_dynamical_correlator_window's
    own comment in chain_session.h) rather than pyitensor's "fix at the
    source" `eshift`-in-matvec approach, since ITensorTDVP's own boundary-
    tensor `tdvp()` overload is a vendored black box that cannot be
    wrapped the same way -- both are mathematically exact (the baseline
    provably factors as a uniform additive-to-identity scalar), verified
    independently on the pyitensor side by a dense-matrix exact-evolution
    cross-check. Two independent (unseeded) `_build_fast(3)` runs of the
    *same* physical model should now give closely-agreeing S(x,t), not
    just a coincidentally-agreeing e0."""
    n_window, dt, nt = 4, 0.1, 5
    cv_a = _build_fast(3)
    cv_b = _build_fast(3)
    assert cv_a.e0 == pytest.approx(cv_b.e0, abs=0.05)

    ts_a, xs_a, S_a = cv_a._session3.td_dynamical_correlator_window(
        n_window, "Sz", "Sz", dt, nt, [0, 1], 30, 1e-10, 50, True, 0)
    ts_b, xs_b, S_b = cv_b._session3.td_dynamical_correlator_window(
        n_window, "Sz", "Sz", dt, nt, [0, 1], 30, 1e-10, 50, True, 0)
    S_a, S_b = np.array(S_a), np.array(S_b)
    # Loose tolerance (independent ground states/gauges, loosely converged
    # -- see _build_fast's own docstring): the point is that these are now
    # close at all, not exact agreement. Pre-fix, this max|diff| was
    # order-1 (comparable to S itself) for the analogous pyitensor check.
    assert np.max(np.abs(S_a - S_b)) < 0.3


def test_td_dynamical_correlator_public_api_dispatches_to_v3():
    """Infinite_Many_Body_Chain.td_dynamical_correlator(itensor_version=3)
    (the public, documented entry point) returns finite S(k,omega) with
    the expected shape, exercising the full dispatch/FFT-reduction path in
    infinitechain.py (not just the raw C++ binding used by the more
    targeted tests above)."""
    cv = _build_fast(3)
    ks = np.linspace(-np.pi, np.pi, 5)
    ks_out, es, Skw = cv.td_dynamical_correlator(
        "Sz", 0, "Sz", n_window=8, dt=0.05, nt=6, x_values=[-1, 0, 1],
        maxdim=40, cutoff=1e-10, niter=100, ks=ks, delta=0.2)
    assert Skw.shape == (len(ks_out), len(es))
    assert np.all(np.isfinite(Skw))
