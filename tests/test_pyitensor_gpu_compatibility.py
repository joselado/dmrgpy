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


def test_tebd_matches_numpy(numpy_backend_restored):
    def build():
        sc = _heisenberg(6)
        sc.tevol_method = "TEBD"
        wf = sc.Sx[0] * sc.get_gs()
        wf = wf * (1.0 / np.sqrt(abs(wf.dot(wf))))
        _, sz = timedependent.evolve_and_measure(sc, operator=sc.Sz[1],
                                                 nt=20, dt=0.05, wf=wf)
        return sz

    ref, got = _on_both(build)
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
