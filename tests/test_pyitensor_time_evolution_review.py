"""Regressions for the 2026-09-26 review of `itensor_version="python"`
real-time evolution.

The review found the integrators themselves right: two-site and one-site
TDVP reproduce expm(-i*dt*H) to ~1e-11 at full bond dimension on spin-1/2
(long range, Dzyaloshinskii-Moriya), spin-1, Jordan-Wigner fermions and
bosons, at real and complex dt, and below full bond dimension they track
ED exactly as closely as v3 does. What it found wrong sat around them:

1. `tdvp._lanczos_expm_multiply`, the Krylov exponentiator every TDVP
   route goes through, stopped on ||v0||*beta_k*|c_k| < 1e-10: no time
   factor, and absolute in the vector. So it carried the units of H and
   the norm of the state (the "absolute Krylov error goal" the 2026-09-25
   record left open), and when its `niter` budget ran out it returned the
   unconverged vector silently.
2. `Chain.evolve_and_measure_tdvp`/`_tebd` never undid the norm the SVD
   truncation removes, where `quench_tdvp` and the v3/julia_live loops do,
   so every <psi(t)|O|psi(t)> carried the weight discarded so far (v3's
   `evolve_and_measure_tebd` had the same defect and was fixed with it).
3. `gse.global_subspace_expand` weighed each Krylov vector H^k|phi> by
   its own norm where v3's `addBasis` normalizes it, so the basis the
   expansion chose depended on the units of H.
"""

import numpy as np
import pytest

from dmrgpy import cppext, spinchain, timedependent
from dmrgpy.pyitensor import tdvp


@pytest.fixture
def krylov_default_restored():
    """Process-wide state, including the device path's k-hint, which a
    forced deferred run leaves wherever its last call stopped: left at a
    large k it makes the next forced run speculate far past its own stop,
    which test_tdz_gpu_batching's matvec-overhead bound measures."""
    hint = tdvp._KRYLOV_K_HINT[0]
    yield
    tdvp.set_krylov_defer_sync(None)
    tdvp._KRYLOV_K_HINT[0] = hint


def _random_hermitian(dim, seed):
    """A dense Hermitian matrix of unit spectral width, and a unit vector."""
    rng = np.random.default_rng(seed)
    X = rng.standard_normal((dim, dim)) + 1j * rng.standard_normal((dim, dim))
    A = (X + X.conj().T) / 2
    w = np.linalg.eigvalsh(A)
    A = A / (w[-1] - w[0])
    v = rng.standard_normal(dim) + 1j * rng.standard_normal(dim)
    return A, v / np.linalg.norm(v)


def _exact(A, v, coeff):
    w, U = np.linalg.eigh(A)
    return U @ (np.exp(coeff * w) * (U.conj().T @ v))


# ------------------------------------------------ 1. Krylov exponentiator

@pytest.mark.parametrize("s", [1.0, 1e-4, 1e-8, 1e-12])
def test_krylov_exponentiator_does_not_depend_on_the_units_of_h(s):
    """s*A evolved for 2/s is one dimensionless problem at every s. With
    the old estimate the error was 1.3e-11 at s=1, 2.6e-3 at s=1e-8 and
    0.49 at s=1e-11 (the Krylov space shrank to one vector)."""
    A, v0 = _random_hermitian(400, seed=0)
    ref = _exact(A, v0, -2j)
    out = tdvp._lanczos_expm_multiply(lambda x: (s * A) @ x, v0, -2j / s,
                                      niter=50)
    assert np.linalg.norm(out - ref) < 1e-10


@pytest.mark.parametrize("c", [1e3, 1.0, 1e-8])
def test_krylov_exponentiator_error_is_relative_to_the_vector(c):
    """A state of norm 1e-8 used to come back 2.6e-3 off, relatively: the
    error goal was absolute in the vector."""
    A, v0 = _random_hermitian(400, seed=1)
    ref = _exact(A, v0, -2j)
    out = tdvp._lanczos_expm_multiply(lambda x: A @ x, c * v0, -2j, niter=50)
    assert np.linalg.norm(out / c - ref) < 1e-10


@pytest.mark.parametrize("defer", [False, True])
def test_krylov_exponentiator_substeps_when_the_budget_runs_out(
        defer, krylov_default_restored):
    """|coeff| times the spectral width = 100 needs more than niter=50
    Lanczos vectors. The unconverged vector used to be returned as it was,
    0.31 off; the step is now split. The deferred (device) path shares
    the sub-stepping and must pick the same sub-steps."""
    A, v0 = _random_hermitian(400, seed=2)
    ref = _exact(A, v0, -100j)
    tdvp.set_krylov_defer_sync(defer)
    out = tdvp._lanczos_expm_multiply(lambda x: A @ x, v0, -100j, niter=50)
    assert np.linalg.norm(out - ref) < 1e-9
    tdvp.set_krylov_defer_sync(False)
    host = tdvp._lanczos_expm_multiply(lambda x: A @ x, v0, -100j, niter=50)
    assert np.max(np.abs(out - host)) < 1e-13


def test_krylov_exponentiator_imaginary_time_and_zero_modes():
    """Real coeff (METTS's imaginary time), where the first column of
    exp(coeff*T) is not bounded by one; and an exact zero mode, which the
    relative exhaustion test must stop on rather than divide by zero."""
    A, v0 = _random_hermitian(400, seed=3)
    for coeff in (-0.1, -20.0):
        ref = _exact(A, v0, coeff)
        out = tdvp._lanczos_expm_multiply(lambda x: A @ x, v0, coeff, niter=50)
        assert np.linalg.norm(out - ref) / np.linalg.norm(ref) < 1e-10
    out = tdvp._lanczos_expm_multiply(lambda x: 0 * x, v0, -1j)
    assert np.all(np.isfinite(out))
    assert np.linalg.norm(out - v0) < 1e-14


# -------------------------------------- 2. and 3., through the public API

def _neel_quench(version, n, s=1.0, method="TDVP", maxm=64):
    """A chain holding s*H (Heisenberg), and the Neel state as its start:
    the ground state of a staggered field, which is exactly a product
    state, so no solver error enters at any s."""
    sc = spinchain.Spin_Chain(["S=1/2"] * n, itensor_version=version)
    h0 = sum((-1) ** i * sc.Sz[i] for i in range(n))
    h1 = sum(s * (sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1]
                  + sc.Sz[i] * sc.Sz[i + 1]) for i in range(n - 1))
    sc.maxm, sc.cutoff, sc.tevol_method = maxm, 1e-12, method
    sc.set_hamiltonian(h0)
    wf = sc.get_gs()
    sc.set_hamiltonian(h1)
    return sc, h1, wf


@pytest.mark.parametrize("method", ["TDVP", "TDVP_GSE"])
def test_real_time_evolution_does_not_depend_on_the_units_of_h(method):
    """The same quench written in units 1e-12 times smaller (dt 1e12 times
    longer) is the same trajectory. On TDVP the Krylov exponentiator broke
    this (0.48 off ED at s=1e-12, 2.5e-4 at 1e-8); on TDVP_GSE the
    expansion did as well (1.0e-4 at s=1, 1.3e-3 at every s<=1e-4)."""
    op = None
    out = {}
    for s in (1.0, 1e-12):
        sc, h1, wf = _neel_quench("python", 8, s=s, method=method)
        op = sc.Sz[0] + sc.Sx[3] * sc.Sy[4]
        _ts, out[s] = timedependent.evolve_and_measure(
                sc, operator=op, nt=21, dt=0.1 / s, wf=wf)
    assert np.max(np.abs(out[1e-12] - out[1.0])) < 1e-8


def test_td_correlator_of_small_operators_is_their_scaled_correlator():
    """C[eps*A, eps*B] = eps^2 C[A, B]. The record's own discriminant: on
    "python" it was 8.76e-2 relative at eps=1e-8 and 2.42e-1 at 1e-9."""
    n = 6
    es = np.linspace(-0.5, 4, 40)

    def spectrum(eps):
        sc = spinchain.Spin_Chain(["S=1/2"] * n, itensor_version="python")
        sc.set_hamiltonian(sum(sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1]
                               + sc.Sz[i] * sc.Sz[i + 1] for i in range(n - 1)))
        sc.maxm = 32
        return sc.get_dynamical_correlator(
                submode="TD", name=(eps * sc.Sz[0], eps * sc.Sz[3]),
                es=es, delta=0.2)[1] / eps ** 2

    ref = spectrum(1.0)
    got = spectrum(1e-9)
    assert np.max(np.abs(got - ref)) / np.max(np.abs(ref)) < 1e-8


@pytest.mark.parametrize("method", ["TDVP", "TEBD"])
def test_evolve_and_measure_restores_the_norm_truncation_removes(method):
    """At maxm=8 a 12-site Neel quench truncates at every step. The norm
    used to drift with it on "python" (<psi|psi> 0.982 by t=5 under TDVP,
    0.970 under TEBD, and <H> 4.9e-2 off its conserved value where v3's
    TDVP drifted 1.8e-4), because nothing put it back."""
    sc, h1, wf = _neel_quench("python", 12, method=method, maxm=8)
    _ts, norms = timedependent.evolve_and_measure(
            sc, operator=sc.Sz[0] * 0 + 1, nt=101, dt=0.05, wf=wf)
    assert np.max(np.abs(norms - 1.0)) < 1e-10
    _ts, es = timedependent.evolve_and_measure(
            sc, operator=h1, nt=101, dt=0.05, wf=wf)
    # what is left is the truncation's own drift of <H>, 1.8e-4 for TDVP
    # and 1.8e-3 for TEBD here
    assert np.max(np.abs(es - es[0])) < 5e-3


@pytest.mark.skipif(not cppext.available(3), reason="needs the v3 extension")
@pytest.mark.parametrize("method", ["TDVP", "TEBD"])
def test_evolve_and_measure_truncated_trajectory_matches_v3(method):
    """Same algorithm, same truncation, same start: python and v3 give the
    same truncated trajectory once both restore the norm (v3's TEBD loop
    did not either, and was fixed with it)."""
    out = {}
    for version in ("python", 3):
        sc, h1, wf = _neel_quench(version, 12, method=method, maxm=8)
        _ts, out[version] = timedependent.evolve_and_measure(
                sc, operator=h1, nt=101, dt=0.05, wf=wf)
    assert np.max(np.abs(out["python"] - out[3])) < 1e-6
