"""The two device-oriented optimizations on the complex-time-evolution
(submode="TDZ") dynamical correlator, checked where they can be checked
without a GPU.

Both changes exist to remove *array dispatches and synchronizations*, not
to change arithmetic (docs/pyitensor_gpu_port_plan.md Sec. 9), so what
needs guarding is precisely that they change no number:

* `tdvp.set_krylov_defer_sync` -- the Krylov exponentiator keeps
  alpha/beta on the device and speculates past its own stopping point,
  then rolls back to the exact k the per-iteration host loop would have
  chosen. The rollback is the risky part and it is pure Python logic, so
  it is tested here by *forcing* the deferred path on NumPy, where the
  unmodified host loop is available as an exact reference. A GPU is not
  needed to catch a wrong k.
* `mpsalgebra.BatchedBras` -- the n_max+1 phi^(n) overlaps as one batched
  sweep instead of n_max+1 separate ones, with the bras zero-padded to a
  common per-bond width.

The device-*residency* half of these claims (which agreement cannot
catch, since a host round trip returns the same numbers) lives in
tests/test_pyitensor_gpu_backend.py, which needs JAX.
"""
import numpy as np
import pytest

from dmrgpy import spinchain
from dmrgpy import tdz as tdzmod
from dmrgpy.pyitensor import tdvp
from dmrgpy.pyitensor.mpsalgebra import BatchedBras, inner


@pytest.fixture
def krylov_default_restored():
    """Process-wide state, so put it back whatever the test does."""
    yield
    tdvp.set_krylov_defer_sync(None)
    tdzmod.set_batched_bras(True)


def _heisenberg(n=8, maxm=24, seed=3):
    np.random.seed(seed)
    sc = spinchain.Spin_Chain(["S=1/2" for _ in range(n)],
                              itensor_version="python")
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1]
        h = h + sc.Sy[i] * sc.Sy[i + 1]
        h = h + sc.Sz[i] * sc.Sz[i + 1]
    sc.set_hamiltonian(h)
    sc.maxm = maxm
    sc.nsweeps = 3
    return sc


def test_deferred_krylov_sync_picks_the_same_stopping_point(krylov_default_restored):
    """The speculative path must stop at the *same* Krylov dimension as
    the per-iteration loop, not merely at a converged one.

    Stopping later would still be accurate (more Krylov vectors is a
    better approximation) and would silently break the cross-backend
    agreement the port promises at 1e-14, so "close enough" is not the
    property to assert here. A dense Hermitian operator is enough: the
    stopping rule sees only alpha/beta.
    """
    rng = np.random.default_rng(4)
    dim = 48
    A = rng.standard_normal((dim, dim)) + 1j * rng.standard_normal((dim, dim))
    A = A + A.conj().T
    # Normalized to unit spectral radius so that a 30-vector Krylov space
    # (the default niter) genuinely converges for every coeff below --
    # otherwise the exact cross-check would be testing the subspace size,
    # not the stopping rule.
    A = A / np.max(np.abs(np.linalg.eigvalsh(A)))
    v0 = rng.standard_normal(dim) + 1j * rng.standard_normal(dim)

    seen = {}

    def matvec_counting(tag):
        def mv(v):
            seen[tag] = seen.get(tag, 0) + 1
            return A @ v
        return mv

    # coeff sweeps a range of step sizes, so the loop converges at
    # several different k rather than always at the same one.
    for coeff in (-0.02j, -0.1j, -0.5j, -2.0j):
        tdvp.set_krylov_defer_sync(False)
        ref = tdvp._lanczos_expm_multiply(matvec_counting("host"), v0, coeff)
        tdvp.set_krylov_defer_sync(True)
        got = tdvp._lanczos_expm_multiply(matvec_counting("dev"), v0, coeff)
        assert np.max(np.abs(got - ref)) < 1e-13, "coeff=%r" % (coeff,)

        evals, evecs = np.linalg.eigh(A)
        exact = evecs @ (np.exp(coeff * evals) * (evecs.conj().T @ v0))
        assert np.max(np.abs(ref - exact)) < 1e-9

    # Speculation costs matvecs; it must not cost an unbounded number of
    # them. _KRYLOV_CHECK_BLOCK-1 per call is the worst case, and the
    # k-hint should keep the realized overhead well under that.
    overhead = seen["dev"] - seen["host"]
    assert 0 <= overhead <= 4 * tdvp._KRYLOV_CHECK_BLOCK


def test_deferred_krylov_sync_survives_lanczos_exhaustion(krylov_default_restored):
    """Speculating past a beta of ~0 must not put nan into the basis.

    A rank-2 operator exhausts its Krylov sequence after two iterations,
    where the host loop returns immediately and the deferred path keeps
    going to its checkpoint -- dividing by that beta is exactly the
    inf/nan trap _lanczos_expm_device clamps against.
    """
    rng = np.random.default_rng(9)
    dim = 24
    u = rng.standard_normal(dim) + 1j * rng.standard_normal(dim)
    w = rng.standard_normal(dim) + 1j * rng.standard_normal(dim)
    u /= np.linalg.norm(u)
    w = w - (u.conj() @ w) * u
    w /= np.linalg.norm(w)
    A = 2.0 * np.outer(u, u.conj()) - 1.5 * np.outer(w, w.conj())
    v0 = u + 0.5 * w

    tdvp.set_krylov_defer_sync(False)
    ref = tdvp._lanczos_expm_multiply(lambda v: A @ v, v0, -0.3j)
    tdvp.set_krylov_defer_sync(True)
    got = tdvp._lanczos_expm_multiply(lambda v: A @ v, v0, -0.3j)
    assert np.all(np.isfinite(got))
    assert np.max(np.abs(got - ref)) < 1e-13


def test_batched_bras_reproduce_inner(krylov_default_restored):
    """One batched sweep against n_max+1 separate ones, on bras with
    genuinely different bond dimensions -- the zero-padding to a common
    per-bond width is what has to be exact."""
    sc = _heisenberg()
    n = sc.ns
    sc.get_gs()
    Hop = sc.toMPO(sc.hamiltonian - sc.e0)
    bras = [sc.toMPO(sc.Sz[n // 2]) * sc.wf0]
    for _ in range(4):
        bras.append(Hop * bras[-1])
    ket = sc.toMPO(sc.Sx[2]) * sc.wf0

    ref = np.array([inner(b.cpp_handle, ket.cpp_handle) for b in bras])
    got = BatchedBras([b.cpp_handle for b in bras]).overlaps(ket.cpp_handle)

    assert got.shape == ref.shape
    scale = max(1.0, float(np.max(np.abs(ref))))
    assert np.max(np.abs(got - ref)) < 1e-12 * scale
    # The bras really do differ in bond dimension, or the padding above is
    # untested rather than exact.
    def bond_widths(chain):
        return tuple(chain.A(i).inds[-1].dim for i in range(1, chain.length()))
    assert len({bond_widths(b.cpp_handle) for b in bras}) > 1


def test_tdz_spectrum_is_unchanged_by_either_optimization(krylov_default_restored):
    """End to end: the same TDZ spectrum with each switch either way.

    This is the claim that actually matters to a user -- the optimizations
    are a scheduling change, so submode="TDZ" must return the same
    correlator however the overlaps are batched and whenever the Krylov
    coefficients are read back."""
    es = np.linspace(-0.5, 3.0, 40)

    def spectrum(defer, batched):
        tdvp.set_krylov_defer_sync(defer)
        tdzmod.set_batched_bras(batched)
        sc = _heisenberg()
        n = sc.ns
        return sc.get_dynamical_correlator(
                mode="DMRG", submode="TDZ",
                name=(sc.Sz[n // 2], sc.Sz[n // 2]),
                dt=0.1, nt=25, es=es)[1]

    ref = spectrum(False, False)
    for defer, batched in ((True, False), (False, True), (True, True)):
        got = spectrum(defer, batched)
        assert np.max(np.abs(got - ref)) < 1e-9, (defer, batched)
