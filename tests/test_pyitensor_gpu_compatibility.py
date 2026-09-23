"""Calculations that had never run on the JAX array backend before
2026-09-16: do they run there at all, and do they give the NumPy answer?

tests/test_pyitensor_gpu_backend.py covers the paths the GPU port was
built around (ground state, static/KPM/TDZ correlators, TDVP). This file
covers the rest of the pure-Python engine, found by running every one of
them once on a GTX 1060 against a NumPy reference (see
docs/gpu_cpu_performance.md's consumer-GPU section).

That sweep found one module that did not run on a device at all:
**iDMRG** (pyitensor/idmrg.py, which had never been ported to
pyitensor/backend.py). It failed twice, in the two ways the backend
module's docstring warns about:

* `np.take` on tensor data. NumPy forwards a free function on a device
  array to `jnp.take`, but with NumPy's default `mode="raise"`, which JAX
  does not implement -- so `gs_energy()` raised before its first growth
  step. Now `bk.xp().take`.
* an in-place `arr[idx] -= ...` in `_subtract_energy_baseline`. JAX arrays
  are immutable. Now `bk.setblock`, which is the same in-place write on
  NumPy.

With both fixed, iDMRG agrees with NumPy to ~1e-12 in energy and ~4e-10
in `vev`/`correlator`/`local_excitation_gap` (an iterative fixed-point
solve, so roundoff from a different BLAS carries through). The other
calculations here ran unchanged and agree to 1e-11..1e-16; they are
pinned so that stays true.

The same sweep run again on 2026-09-22, for the last two unported
modules (docs/pyitensor_gpu_port_plan.md Sec. 9 item 6), found the other
shape of the same problem in **GSE** (pyitensor/gse.py): it ran on a
device and gave the right answer, and rebuilt one tensor per bond on the
host while doing it, because `np.concatenate` on device inputs returns a
host array rather than raising. Agreement cannot see that, so the two
residency tests at the end of this file assert on the array type
instead. TEBD (pyitensor/tebd.py) was measured the same way and needed
nothing: its NumPy calls build the bond gates once at setup and the
evolution itself never leaves the device.

Skipped when JAX is not installed; not skipped on a CPU-only JAX, for the
reason tests/test_pyitensor_gpu_backend.py gives (same code path, same
immutability, just slower). Sizes are tiny for the same reason too.
"""

import numpy as np
import pytest

from dmrgpy import fermionchain, infinitechain, spinchain, timedependent

jax = pytest.importorskip("jax", reason="the JAX backend needs jax installed")

from dmrgpy.pyitensor import backend as bk       # noqa: E402


@pytest.fixture
def numpy_backend_restored():
    """The backend is process-wide state: put it back whatever the test
    does, or the rest of the suite silently runs on JAX."""
    yield
    bk.set_backend("numpy")
    bk.set_pad_bonds(None)
    bk.set_jit("auto")


def _on_both(build):
    """`build()` on NumPy, then on JAX, from the same random seed."""
    np.random.seed(1234)
    ref = build()
    bk.set_backend("jax")
    np.random.seed(1234)
    got = build()
    return np.asarray(ref, dtype=complex), np.asarray(got, dtype=complex)


def _heisenberg(n, maxm=16, nsweeps=6):
    sc = spinchain.Spin_Chain(["S=1/2"] * n, itensor_version="python")
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] \
            + sc.Sz[i] * sc.Sz[i + 1]
    sc.set_hamiltonian(h)
    sc.maxm = maxm
    sc.nsweeps = nsweeps
    return sc


def _interacting_fermions(n, maxm=16, nsweeps=6):
    fc = fermionchain.Fermionic_Chain(n, itensor_version="python")
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdag[i] * fc.C[i + 1]
    h = h + h.get_dagger()
    for i in range(n - 1):
        h = h + (fc.N[i] - 0.5) * (fc.N[i + 1] - 0.5)
    for i in range(n):
        h = h + 0.3 * ((-1) ** i) * fc.N[i]
    fc.set_hamiltonian(h)
    fc.maxm = maxm
    fc.nsweeps = nsweeps
    return fc


def _dimerized_infinite_chain(gs_method, maxm):
    """Gapped (dimerized, staggered field), so both solvers converge fast
    and every observable below is nonzero -- `<Sz>=0` on a uniform chain
    would pass however wrong the gauge was."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"],
                                           itensor_version="python")
    ic.gs_method = gs_method
    h = (ic.SxC[0] * ic.SxC[1] + ic.SyC[0] * ic.SyC[1] + ic.SzC[0] * ic.SzC[1]
         + 0.4 * (ic.SxC[1] * ic.SxR[0] + ic.SyC[1] * ic.SyR[0]
                  + ic.SzC[1] * ic.SzR[0])
         + 0.2 * ic.SzC[0] - 0.2 * ic.SzC[1])
    ic.maxm = maxm
    ic.maxiter = 60
    ic.set_hamiltonian(h)
    return ic


def _infinite_observables(ic):
    return [ic.gs_energy(), ic.vev("Sz", 0), ic.vev("Sz", 1),
            ic.correlator("Sz", 0, "Sz", 1), ic.correlator("Sx", 0, "Sx", 3)]


def test_idmrg_runs_on_jax_and_matches_numpy(numpy_backend_restored):
    """The regression for both idmrg.py fixes: before them this raised
    inside gs_energy() on any device."""
    def build():
        ic = _dimerized_infinite_chain("idmrg", maxm=12)
        return _infinite_observables(ic) + [ic.local_excitation_gap()]

    ref, got = _on_both(build)
    assert abs(ref[0].imag) < 1e-12 and ref[1].real != pytest.approx(0.0)
    assert got[0] == pytest.approx(ref[0], abs=1e-10)          # energy
    assert np.max(np.abs(got[1:5] - ref[1:5])) < 1e-8           # vev, corr
    assert got[5] == pytest.approx(ref[5], abs=1e-7)            # local gap


def test_vumps_matches_numpy(numpy_backend_restored):
    ref, got = _on_both(
        lambda: _infinite_observables(_dimerized_infinite_chain("vumps", maxm=8)))
    assert np.max(np.abs(got - ref)) < 1e-9


def test_conserved_sector_ground_state_matches_numpy(numpy_backend_restored):
    """Sector mode on this backend is dense storage plus a charge penalty on
    the variational solve (pyitensor/sector.py) -- a separate code path
    from the plain ground state."""
    def build():
        fc = _interacting_fermions(6)
        fc.set_conserved_sector(Nf=3)
        return [fc.gs_energy()]

    ref, got = _on_both(build)
    assert got[0] == pytest.approx(ref[0], abs=1e-9)


def _quench_sz(n, tevol_method, nt=20, dt=0.05, **chain_attrs):
    """<Sz_0>(t) after preparing the ground state of a staggered field
    (plus a little XY, so it is not a product state) and quenching to
    Heisenberg.

    The observable has to be chosen with some care here, and the obvious
    choice is empty: evolving Sx[0]|GS> of a uniform Heisenberg chain and
    measuring Sz[1] gives identically zero, because that state is
    invariant under a global pi rotation about x, which flips Sz and
    commutes with H. Measured, the whole trajectory sits at ~1e-13, so
    a cross-backend comparison of it compares roundoff to roundoff and
    an integrator that had corrupted the state would still pass. The
    quench below carries <Sz_0> ~ -0.49 instead; measured on a 30-step
    version of it, the two array backends agree at 1.2e-14 for TEBD and
    7.8e-15 for TDVP_GSE, which is what sets the 1e-9 the callers assert.
    """
    sc = spinchain.Spin_Chain([2] * n, itensor_version="python")
    sc.tevol_method = tevol_method
    sc.maxm = 20
    sc.nsweeps = 8
    for name, value in chain_attrs.items():
        setattr(sc, name, value)

    h0 = 0
    for i in range(n):
        h0 = h0 + (-1) ** i * sc.Sz[i]
    for i in range(n - 1):
        h0 = h0 + 0.3 * (sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1])
    h1 = 0
    for i in range(n - 1):
        h1 = h1 + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] \
            + sc.Sz[i] * sc.Sz[i + 1]

    sc.set_hamiltonian(h0)
    wf = sc.get_gs()
    sc.set_hamiltonian(h1)
    _, sz = timedependent.evolve_and_measure(sc, operator=sc.Sz[0],
                                             nt=nt, dt=dt, wf=wf)
    return sz


def test_tebd_matches_numpy(numpy_backend_restored):
    ref, got = _on_both(lambda: _quench_sz(6, "TEBD"))
    assert np.max(np.abs(ref.real)) > 0.4, "the observable is trivially zero"
    assert np.max(np.abs(got - ref)) < 1e-9


def test_batched_four_point_tensor_matches_numpy(numpy_backend_restored):
    """ctmode="batched" is the kernel docs/gpu_cpu_performance.md singles
    out as winning on a device at small bond dimension -- worth pinning
    that it is also right there."""
    def build():
        wf = _interacting_fermions(4).get_gs()
        return np.ravel(wf.get_four_correlation_tensor(ctmode="batched"))

    ref, got = _on_both(build)
    assert np.max(np.abs(got - ref)) < 1e-10


def test_non_hermitian_dmrg_matches_numpy(numpy_backend_restored):
    def build():
        fc = fermionchain.Fermionic_Chain(6, itensor_version="python")
        h = 0
        for i in range(5):
            h = h + 1.2 * fc.Cdag[i] * fc.C[i + 1] + 0.8 * fc.Cdag[i + 1] * fc.C[i]
            h = h + (fc.N[i] - 0.5) * (fc.N[i + 1] - 0.5)
        fc.set_hamiltonian(h)
        fc.maxm = 20
        fc.nsweeps = 8
        e, _psil, _psir = fc.nhdmrg()
        return [e]

    ref, got = _on_both(build)
    assert got[0] == pytest.approx(ref[0], abs=1e-9)


def _itensor_source_types(run):
    """The array type every ITensor built during `run()` was constructed
    from, as a Counter.

    This is the only way to see the class of host transfer
    docs/pyitensor_gpu_port_plan.md Sec. 5 warns about. A free NumPy
    function applied to a device array, `np.concatenate` being the one
    that mattered here, returns a *host* array with no error and no
    exception, and the next ITensor built from it silently converts back:
    the numbers are identical and only the time changes. A `backend.
    to_host` counter cannot see it either, since none of those calls goes
    through `to_host`. The array type does.
    """
    from collections import Counter

    from dmrgpy.pyitensor.tensor import ITensor

    seen = []
    real_init = ITensor.__init__

    def recording_init(self, inds, array=None):
        real_init(self, inds, array)
        if array is not None:
            # isinstance, not type(...).__name__: jaxlib's concrete array
            # class has already moved module once and its name is not API,
            # while jax.Array is.
            seen.append("host" if isinstance(array, np.ndarray)
                        else "device" if isinstance(array, jax.Array)
                        else type(array).__name__)

    ITensor.__init__ = recording_init
    try:
        run()
    finally:
        ITensor.__init__ = real_init
    return Counter(seen)


def _gse_sweep_on(sites, terms, chi=8):
    """One global_subspace_expand() call on a random MPS, which is the
    whole of pyitensor/gse.py's array work."""
    from dmrgpy.pyitensor.autompo import AutoMPO
    from dmrgpy.pyitensor.gse import global_subspace_expand
    from dmrgpy.pyitensor.mpobuilder import to_mpo
    from dmrgpy.pyitensor.mpsalgebra import randomMPS

    H = to_mpo(AutoMPO.from_terms(sites, terms))
    psi = randomMPS(sites, chi)
    psi.position(1)
    return lambda: global_subspace_expand(H, psi, 3, 1e-8, maxdim=20,
                                          bond_maxdim=20)


def _heisenberg_terms(n):
    return [(1.0, [(op, i), (op, i + 1)])
            for i in range(1, n) for op in ("Sx", "Sy", "Sz")]


def test_tdvp_gse_matches_numpy(numpy_backend_restored):
    """tevol_method="TDVP_GSE" is one-site TDVP plus pyitensor/gse.py's
    Krylov basis enrichment for the leading tdvp_gse_sweeps steps, so it
    is the only route that exercises gse.py at all. tdvp_gse_sweeps is 5
    of the 20 steps rather than the default 3, so the expansion runs on a
    state that has already evolved."""
    ref, got = _on_both(lambda: _quench_sz(
        6, "TDVP_GSE", tdvp_gse_sweeps=5, tdvp_gse_krylov_order=3,
        tdvp_gse_cutoff=1e-8))
    assert np.max(np.abs(ref.real)) > 0.4, "the observable is trivially zero"
    assert np.max(np.abs(got - ref)) < 1e-9


def test_gse_keeps_every_tensor_on_the_device(numpy_backend_restored):
    """The regression for gse.py's port. `np.concatenate` of V1's rows
    with U2's new directions returned a host array from device inputs, so
    the enlarged tensor res.A(b) was rebuilt on the host once per bond:
    measured at 5 host-sourced tensors out of the 229 a 6-site expansion
    builds, exactly one per bond. Agreement cannot see that, which is why
    this asserts on the array type instead."""
    from dmrgpy.pyitensor.sites import SiteX

    n = 6
    bk.set_backend("jax")
    sites = SiteX([2] * n)   # 2 = SpinHalfSite (siteset.TYPE_CODE_TO_SITE)
    counts = _itensor_source_types(_gse_sweep_on(sites, _heisenberg_terms(n)))
    assert sum(counts.values()) > 100, counts
    assert set(counts) == {"device"}, counts


def test_tebd_keeps_every_tensor_on_the_device(numpy_backend_restored):
    """tebd.py was never ported and did not need to be: its NumPy calls
    build the bond Hamiltonians and exponentiate them with scipy, once at
    setup, and the gate crosses to the device at ITensor.__init__ like
    any other array. What matters is that the *evolution* adds no further
    crossing, so this pins a step() rather than the setup."""
    from dmrgpy.pyitensor.mpsalgebra import randomMPS
    from dmrgpy.pyitensor.sites import SiteX
    from dmrgpy.pyitensor.tebd import TEBDEvolver

    n = 6
    bk.set_backend("jax")
    sites = SiteX([2] * n)
    psi = randomMPS(sites, 8)
    psi.position(1)
    evolver = TEBDEvolver(sites, _heisenberg_terms(n), 0.05, 1e-10, 20)
    gates = list(evolver._gates_half.values()) + list(evolver._gates_full.values())
    assert gates and all(isinstance(g.array, jax.Array) for g in gates)

    counts = _itensor_source_types(lambda: evolver.step(psi))
    assert sum(counts.values()) > 100, counts
    assert set(counts) == {"device"}, counts
