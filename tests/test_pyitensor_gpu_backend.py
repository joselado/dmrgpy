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
    bk.set_jit("auto")


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


def _quench_trajectory(n=6, nt=10, dt=0.05, maxm=20, seed=5):
    """Real-time TDVP after a quench: prepare the ground state of a
    staggered field, evolve under Heisenberg, measure <Sz_0>(t)."""
    from dmrgpy import timedependent
    np.random.seed(seed)
    sc = spinchain.Spin_Chain(["S=1/2" for _ in range(n)],
                              itensor_version="python")
    h0 = 0
    for i in range(n):
        h0 = h0 + (-1) ** i * sc.Sz[i]
    h1 = 0
    for i in range(n - 1):
        h1 = h1 + sc.Sx[i] * sc.Sx[i + 1]
        h1 = h1 + sc.Sy[i] * sc.Sy[i + 1]
        h1 = h1 + sc.Sz[i] * sc.Sz[i + 1]
    sc.set_hamiltonian(h0)
    sc.maxm = maxm
    sc.nsweeps = 3
    wf = sc.get_gs()
    sc.set_hamiltonian(h1)
    return timedependent.evolve_and_measure(sc, operator=sc.Sz[0],
                                            nt=nt, dt=dt, wf=wf)


def test_time_evolution_matches_numpy(numpy_backend_restored):
    """TDVP on the device. Unlike the ground state, real-time evolution
    grows entanglement with time, so this is the calculation whose chi
    climbs into the range where a device wins by construction -- see
    docs/pyitensor_gpu_port_plan.md Sec. 9.

    A whole trajectory rather than one number: the Krylov propagator runs
    once per bond per time step, so an error in its device/host split
    would show up as drift, which comparing only the final point could
    hide."""
    _ts, sz_np = _quench_trajectory()
    bk.set_backend("jax")
    _ts, sz_jax = _quench_trajectory()
    assert np.max(np.abs(sz_jax - sz_np)) < 1e-10


def test_krylov_propagator_keeps_its_basis_on_the_device(numpy_backend_restored):
    """The point of porting tdvp.py was residency, not agreement: the
    Krylov basis is chi^2 d^2 per vector times up to `niter` vectors, and
    before the port it lived in a preallocated NumPy buffer, so every
    iteration pulled the propagated vector back to the host and pushed it
    out again.

    Agreement alone cannot catch a regression to that -- a host round trip
    returns the same numbers, only slower. So this asserts on the array
    *type* that comes out: a device array can only be produced by a basis
    that stayed on the device."""
    from dmrgpy.pyitensor import tdvp

    bk.set_backend("jax")
    rng = np.random.default_rng(11)
    dim = 32
    A = rng.standard_normal((dim, dim)) + 1j * rng.standard_normal((dim, dim))
    A = A + A.conj().T
    A_dev = bk.asarray(A)
    v0 = bk.asarray(rng.standard_normal(dim) + 1j * rng.standard_normal(dim))

    out = tdvp._lanczos_expm_multiply(lambda v: A_dev @ v, v0, -0.1j)
    assert isinstance(out, jax.Array)

    # ...and it is still the right answer: exp(-0.1i A) v0, done densely.
    evals, evecs = np.linalg.eigh(A)
    exact = evecs @ (np.exp(-0.1j * evals) * (evecs.conj().T @ np.asarray(bk.to_host(v0))))
    assert np.max(np.abs(np.asarray(bk.to_host(out)) - exact)) < 1e-9


def test_jit_is_tied_to_padding_and_changes_no_result(numpy_backend_restored):
    """backend.set_jit's "auto" rule: fuse the hot composites exactly when
    set_pad_bonds has frozen their shapes, since jax.jit traces per shape
    and DMRG mints a new one at every bond dimension. Both knobs are pure
    performance, so the energy must not move under any combination."""
    assert not bk.jit_enabled()          # never on the NumPy backend
    bk.set_backend("jax")
    assert not bk.jit_enabled()          # auto, and nothing is padded yet
    bk.set_pad_bonds(16)
    assert bk.jit_enabled()              # auto + padding -> on
    bk.set_jit(False)
    assert not bk.jit_enabled()          # explicit override wins
    e_eager = float(np.real(_heisenberg().gs_energy(mode="DMRG")))
    bk.set_jit(True)
    assert bk.jit_enabled()
    e_jit = float(np.real(_heisenberg().gs_energy(mode="DMRG")))
    assert e_jit == pytest.approx(e_eager, abs=1e-9)
    comps = bk.compilations()            # None if JAX drops the private
    assert comps is None or comps > 0    # cache-size API this reads
    with pytest.raises(ValueError):
        bk.set_jit("sometimes")
    bk.set_jit("auto")


def test_svd_keeps_its_singular_value_tensor_on_the_device(numpy_backend_restored):
    """svd() runs at every bond of every half-sweep, so its own output is
    the engine's most frequent opportunity to leak an O(chi^2) transfer.

    The diagonal S tensor used to be assembled with np.diag() from the
    host copy of the spectrum -- correct, and a full keep x keep matrix
    pushed back across the bus per factorization. Only the O(chi)
    spectrum needs to be on the host; this asserts the matrix is not."""
    from dmrgpy.pyitensor.index import Index
    from dmrgpy.pyitensor.svd import svd
    from dmrgpy.pyitensor.tensor import ITensor

    bk.set_backend("jax")
    rng = np.random.default_rng(3)
    a, b = Index(6, tags="Link"), Index(8, tags="Link")
    T = ITensor((a, b), bk.asarray(rng.standard_normal((6, 8))
                                   + 1j * rng.standard_normal((6, 8))))
    U, S, V, spectrum = svd(T, [a], cutoff=1e-12)
    assert isinstance(S.array, jax.Array)
    assert isinstance(U.array, jax.Array) and isinstance(V.array, jax.Array)
    # The spectrum itself is host-side bookkeeping, by design.
    assert not isinstance(spectrum.eigs, jax.Array)
    # ...and the split still reconstructs the tensor it came from.
    back = np.asarray(bk.to_host((U * S * V).transpose_to((a, b))))
    assert np.max(np.abs(back - np.asarray(bk.to_host(T.array)))) < 1e-12


def test_batched_bras_stay_on_the_device(numpy_backend_restored):
    """mpsalgebra.BatchedBras exists to give the phi^(n) overlaps of a
    submode="TDZ" run a batch axis, which is only worth anything if the
    prepared bras live where the ket does -- a batch that transfers per
    step would be strictly worse than the loop it replaces.

    Agreement is checked on NumPy in tests/test_tdz_gpu_batching.py; what
    needs a device is residency, so this asserts on the array type."""
    from dmrgpy.pyitensor.mpsalgebra import BatchedBras, inner

    bk.set_backend("jax")
    sc = _heisenberg(n=4, maxm=8, nsweeps=2)
    sc.get_gs()
    Hop = sc.toMPO(sc.hamiltonian - sc.e0)
    bras = [sc.toMPO(sc.Sz[1]) * sc.wf0]
    for _ in range(2):
        bras.append(Hop * bras[-1])
    ket = sc.toMPO(sc.Sx[0]) * sc.wf0

    bundle = BatchedBras([b.cpp_handle for b in bras])
    assert all(isinstance(arr, jax.Array) for arr in bundle._bras)

    got = bundle.overlaps(ket.cpp_handle)
    ref = np.array([inner(b.cpp_handle, ket.cpp_handle) for b in bras])
    assert np.max(np.abs(got - ref)) < 1e-11


def test_tdz_correlator_matches_numpy(numpy_backend_restored):
    """End to end on the device: complex-time evolution (submode="TDZ",
    arXiv:2311.10909) exercises the deferred-synchronization Krylov path
    and the batched phi^(n) overlaps together, on top of TDVP.

    Tiny on purpose -- eager JAX is far slower than NumPy at this size
    (see this module's docstring); the point is that the two backends
    agree, not that either is fast here."""
    es = np.linspace(-0.5, 2.5, 25)

    def spectrum():
        sc = _heisenberg(n=4, maxm=8, nsweeps=2)
        return sc.get_dynamical_correlator(
                mode="DMRG", submode="TDZ", name=(sc.Sz[1], sc.Sz[1]),
                dt=0.1, nt=8, es=es)[1]

    ref = spectrum()
    bk.set_backend("jax")
    got = spectrum()
    assert np.max(np.abs(got - ref)) < 1e-9
