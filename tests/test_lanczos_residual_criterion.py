"""Coverage for `pyitensor/dmrg.py::_lanczos_ground_state`'s opt-in Ritz
RESIDUAL stopping criterion, and for the VUMPS convergence it unblocks.

The default criterion stops when the lowest Ritz *value* stops moving.
For a Hermitian operator the value error is quadratic in the eigenvector
error, so a value tolerance of `t` leaves the eigenvector accurate only to
~sqrt(t). That is the right trade for finite DMRG, which reports energies,
and the wrong one for VUMPS, whose own convergence criterion is a norm
difference between two independently-solved eigenvectors (AC and C) -- it
made `gauge_mismatch` floor at ~1e-6, so `tol=1e-10` was unreachable at any
number of outer iterations and `.converged` was effectively always False
for D>=8.

The tests below pin, in order: that the two criteria really do differ in
eigenvector accuracy by the square root the theory predicts; that the
default path is untouched; and that VUMPS now actually converges on a case
that previously could not.
"""
import numpy as np
import pytest
from scipy.integrate import quad

from dmrgpy import infinitechain
from dmrgpy.pyitensor import vumps
from dmrgpy.pyitensor.dmrg import _lanczos_ground_state


def _hermitian(n, seed):
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    return (A + A.conj().T) / 2


def _residual(H, lam, v):
    return float(np.linalg.norm(H @ v - lam * v))


def test_residual_criterion_reaches_an_eigenvector_the_value_criterion_cannot():
    """Same operator, same Krylov budget, same starting vector: only the
    stopping rule differs. The residual criterion has to deliver an
    eigenvector orders of magnitude better -- that gap is the entire reason
    the parameter exists."""
    n = 48
    H = _hermitian(n, 11)
    v0 = np.random.default_rng(0).standard_normal(n) + 0j
    # a full Krylov budget, so the stopping rule -- not the iteration cap --
    # is what decides where each one lands
    lam_a, va = _lanczos_ground_state(lambda x: H @ x, v0.copy(), niter=n)
    lam_b, vb = _lanczos_ground_state(lambda x: H @ x, v0.copy(), niter=n,
                                      residual_tol=1e-12)
    ra, rb = _residual(H, lam_a, va), _residual(H, lam_b, vb)
    assert rb < 1e-10, "residual criterion did not deliver its own tolerance"
    assert ra > 1e-8, "value criterion unexpectedly gave an exact eigenvector"
    # And the gap is the square root the theory predicts, not an arbitrary
    # constant: a value tolerance of 1e-12 buys an eigenvector good to
    # ~sqrt(1e-12)=1e-6, which is exactly where `ra` lands.
    assert 1e-8 < ra < 1e-4
    # both still find the same eigenVALUE -- the value was never the problem
    exact = float(np.linalg.eigvalsh(H)[0])
    assert lam_a == pytest.approx(exact, abs=1e-9)
    assert lam_b == pytest.approx(exact, abs=1e-12)


def test_default_path_is_unchanged_by_the_new_parameter():
    """`residual_tol=None` must be the pre-existing code path exactly --
    every finite-DMRG caller relies on that."""
    n = 48
    H = _hermitian(n, 11)
    v0 = np.random.default_rng(5).standard_normal(n) + 0j
    lam_a, va = _lanczos_ground_state(lambda x: H @ x, v0.copy(), niter=30)
    lam_b, vb = _lanczos_ground_state(lambda x: H @ x, v0.copy(), niter=30,
                                      residual_tol=None)
    assert lam_a == lam_b
    assert np.array_equal(np.asarray(va), np.asarray(vb))


def test_converged_warm_start_returns_after_one_matvec():
    """The residual check sits before the loop body precisely so an
    already-converged start costs a single matvec -- which is where a late
    VUMPS iteration lives."""
    n = 32
    H = _hermitian(n, 7)
    w, V = np.linalg.eigh(H)
    calls = [0]

    def matvec(x):
        calls[0] += 1
        return H @ x

    lam, _v = _lanczos_ground_state(matvec, V[:, 0].astype(complex),
                                    niter=30, residual_tol=1e-10)
    assert lam == pytest.approx(float(w[0]), abs=1e-12)
    assert calls[0] == 1


def _tfim(g):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
    ic.set_hamiltonian(-4.0 * ic.SxC[0] * ic.SxR[0] - 2.0 * g * ic.SzC[0])
    return ic


def _tfim_exact_energy_density(g):
    val, _ = quad(lambda k: np.sqrt(1 + g ** 2 - 2 * g * np.cos(k)), 0, np.pi)
    return -val / np.pi


@pytest.mark.parametrize("g", [1.5, 1.0])
def test_vumps_actually_converges_at_D8(g):
    """The end-to-end consequence, and the regression this file exists to
    guard: a D=8 TFIM chain reports `.converged` and a gauge mismatch
    inside the requested tolerance. Before the residual criterion this was
    False for every attempt, at both couplings, however many outer
    iterations were spent -- while the energy was already correct, which is
    why an energy-only test could not see it. Both are asserted here."""
    ic = _tfim(g)
    np.random.seed(0)
    res = vumps.vumps_ground_state(ic.site_types, ic._h_intra.op,
                                   ic._h_inter.op, ic.n_uc, 8,
                                   tol=1e-10, maxiter=800, nrestarts=2)
    assert res.converged
    assert res.gauge_mismatch < 1e-10
    assert res.e0 == pytest.approx(_tfim_exact_energy_density(g), abs=1e-4)
