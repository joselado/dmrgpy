"""The pure-Python engine's array backend (pyitensor/backend.py): does
running on JAX instead of NumPy give the same physics?

Skipped entirely when JAX is not installed. It is NOT skipped when JAX has
only a CPU device: the point of these tests is that the *array library*
swap is faithful, and that is worth checking wherever it can run -- a
host-only JAX run exercises exactly the same code path the device does
(same immutability, same dispatch, same host round trips), just slower.

Sizes are deliberately tiny. Eager JAX is 8-14x slower than NumPy at small
bond dimension -- the device only wins above chi ~ 120-160, see
docs/gpu_cpu_performance.md -- and XLA compiles a kernel per distinct
(operation, shape), so a larger chain here would cost minutes for no extra
coverage.

Every test restores the NumPy backend afterwards through the fixture: the
backend is process-wide state, so leaking it would silently move the rest
of the suite onto JAX.
"""

import numpy as np
import pytest

from dmrgpy import spinchain

jax = pytest.importorskip("jax", reason="the JAX backend needs jax installed")

from dmrgpy.pyitensor import backend as bk       # noqa: E402


@pytest.fixture
def numpy_backend_restored():
    """Process-wide state, so put it back whatever the test does."""
    yield
    bk.set_backend("numpy")
    bk.set_pad_bonds(None)


def _heisenberg(n=6, maxm=20, kpmmaxm=20, nsweeps=4, seed=17):
    # Same seed for both backends: DMRG starts from a random MPS, so a
    # cross-backend comparison is only meaningful from the same start.
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
    sc.kpmmaxm = kpmmaxm
    sc.nsweeps = nsweeps
    return sc


def test_backend_switch_reports_itself(numpy_backend_restored):
    """set_backend is explicit and never falls back silently -- a GPU run
    that quietly became a CPU run is a benchmark that lies."""
    assert bk.backend_name() == "numpy"
    assert not bk.is_device()
    bk.set_backend("jax")
    assert bk.backend_name() == "jax"
    assert "jax" in bk.device_info()
    with pytest.raises(ValueError):
        bk.set_backend("tensorflow")


def test_ground_state_energy_matches_numpy(numpy_backend_restored):
    """The variational energy is the sharpest scalar check available: it is
    sensitive to every tensor in the state. Measured agreement on larger
    chains is ~1e-15 to 1e-11; 1e-9 here leaves room for the iterative
    solve without being loose enough to hide a real error."""
    e_np = float(np.real(_heisenberg().gs_energy(mode="DMRG")))
    bk.set_backend("jax")
    e_jax = float(np.real(_heisenberg().gs_energy(mode="DMRG")))
    assert e_jax == pytest.approx(e_np, abs=1e-9)


def test_static_correlator_matches_numpy(numpy_backend_restored):
    """A two-site expectation value, i.e. the measurement path (device ->
    host for the returned scalar) as well as the sweep."""
    sc = _heisenberg()
    v_np = complex(sc.vev(sc.Sz[1] * sc.Sz[2])).real
    bk.set_backend("jax")
    sc = _heisenberg()
    v_jax = complex(sc.vev(sc.Sz[1] * sc.Sz[2])).real
    assert v_jax == pytest.approx(v_np, abs=1e-9)


def test_kpm_correlator_obeys_the_sum_rule_on_jax(numpy_backend_restored):
    """The KPM moment recursion is the heaviest consumer of the engine
    (applyMPO + MPS sum + truncation per moment), so it exercises svd.py's
    device/host split -- the singular values come back to the host for the
    truncation rule while U and V stay put. Checked against the exact
    zeroth-moment sum rule rather than against NumPy pointwise: the
    recursion truncates ~150 times, so tiny differences in truncation
    decisions accumulate (measured ~1e-5 between backends, and ~1e-14
    between two NumPy runs), while mu_0 = <GS|Sz Sz|GS> = 1/4 is exact for
    a spin-1/2 site whatever path got there."""
    es = np.linspace(-8.0, 8.0, 1500)
    bk.set_backend("jax")
    sc = _heisenberg(kpmmaxm=24)
    x, y = sc.get_dynamical_correlator(mode="DMRG", submode="KPM",
                                        name=(sc.Sz[0], sc.Sz[0]),
                                        es=es, delta=0.15)
    weight = np.trapezoid(np.real(y), x)
    assert weight == pytest.approx(0.25, abs=5e-3)


def test_pad_bonds_leaves_the_ground_state_energy_alone(numpy_backend_restored):
    """set_pad_bonds freezes every bond at one dimension so XLA compiles one
    kernel per operation instead of one per bond dimension. It appends
    *zero* singular values after truncation has already chosen what to
    keep, so the represented state is unchanged and the energy must not
    move. (It does perturb later truncation *decisions* in a long
    recursion -- see docs/gpu_cpu_performance.md -- which is why this
    asserts on the variational energy and not on a KPM spectrum.)"""
    e_plain = float(np.real(_heisenberg().gs_energy(mode="DMRG")))
    bk.set_pad_bonds(20)
    e_padded = float(np.real(_heisenberg().gs_energy(mode="DMRG")))
    assert e_padded == pytest.approx(e_plain, abs=1e-9)


def test_setblock_is_functional_for_immutable_arrays(numpy_backend_restored):
    """backend.setblock exists because JAX arrays cannot be written in
    place, and mpsalgebra.py's direct-sum construction needs a block write.
    Callers must use the return value; this pins that contract for both
    backends."""
    for name in ("numpy", "jax"):
        bk.set_backend(name)
        arr = bk.zeros((4, 4))
        block = bk.asarray(np.ones((2, 2)))
        out = bk.setblock(arr, (slice(0, 2), slice(0, 2)), block)
        assert complex(bk.to_host(out).sum()) == pytest.approx(4.0)
        assert bk.scalar(out[0, 0]) == pytest.approx(1.0)


def test_krylov_energy_projection_matches_across_backends(numpy_backend_restored):
    """pyitensor/kpm_energy_truncation.py's _local_krylov_projection is the
    KPM energy-truncation kernel (Holzner et al., PRB 83, 195115, Sec.
    III-B): build a Krylov basis, diagonalize the projected matrix, and
    drop every component whose Ritz value exceeds a threshold.

    Ported, it splits across the device boundary -- the O(dim) Krylov
    vectors and the final reconstruction stay on the device, while the
    k x k projected matrix, its eigendecomposition and the keep/drop mask
    (k <= dK <= 30) run on the host. This tests that split directly with a
    synthetic Hermitian matvec, rather than through a full KPM run: the
    recursion applies this hundreds of times and truncates in between, so
    an end-to-end comparison is dominated by truncation-trajectory
    divergence (measured: two NumPy runs of the same narrow-window setup
    differ by 1.8e-2 to 2.5e-2 pointwise, while NumPy vs JAX differ by
    5.8e-3 -- i.e. the backends agree better than one backend agrees with
    itself, so the end-to-end number cannot resolve a real error here).
    """
    from dmrgpy.pyitensor import kpm_energy_truncation as et

    dim = 24
    rng = np.random.default_rng(5)
    A = rng.standard_normal((dim, dim)) + 1j * rng.standard_normal((dim, dim))
    A = (A + A.conj().T) / (4 * dim)          # Hermitian, spectrum inside ~[-1,1]
    x0_host = rng.standard_normal(dim) + 1j * rng.standard_normal(dim)

    results = {}
    for name in ("numpy", "jax"):
        bk.set_backend(name)
        A_dev, x0 = bk.asarray(A), bk.asarray(x0_host)
        new_x0, dropped = et._local_krylov_projection(
            lambda v: A_dev @ v, x0, 8, 0.02)
        results[name] = (bk.to_host(new_x0), dropped)

    v_np, w_np = results["numpy"]
    v_jax, w_jax = results["jax"]
    assert w_jax == pytest.approx(w_np, rel=1e-8)
    assert np.max(np.abs(v_jax - v_np)) < 1e-10
    # and the projection actually did something: some weight was removed,
    # so this is not a vacuous comparison of two identity operations
    assert w_np > 0.0
