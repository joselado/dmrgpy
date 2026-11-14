"""TASK (pyitensor/task.py) on the JAX array backend: the Krylov basis
lives on the backend and the spectrum is the NumPy one.

Skipped when JAX is not installed, not when it has only a CPU device, for
the reason test_pyitensor_gpu_backend.py gives: the array library swap is
what is being checked, and a host-only JAX run takes the same code path.
Skipped as well when JAX cannot use the device it finds, for instance a
GPU whose memory another process holds: that is the state of the machine,
not of the code. The probe asks for the devices and allocates one array,
since an out-of-memory can surface at the first allocation rather than at
jax.devices().
"""

import numpy as np
import pytest

from test_kondo_spectrum_task import EV, KW, impurity_chain

jax = pytest.importorskip("jax", reason="the JAX backend needs jax installed")

from dmrgpy.pyitensor import backend as bk  # noqa: E402
from dmrgpy.pyitensor import task  # noqa: E402

# what JAX raises when a device cannot be initialized or allocated on;
# JaxRuntimeError is not a RuntimeError subclass in every release
_DEVICE_ERRORS = tuple({RuntimeError,
                        getattr(getattr(jax, "errors", None),
                                "JaxRuntimeError", RuntimeError)})


def _device_error():
    """None if JAX can allocate on its default device, else the error."""
    try:
        jax.devices()
        (jax.numpy.zeros(16).reshape(4, 4) + 1.0).block_until_ready()
    except _DEVICE_ERRORS as err:  # a GPU held by another process
        return err
    return None


@pytest.fixture(scope="module")
def jax_device():
    err = _device_error()
    if err is not None:
        pytest.skip("JAX could not use its device: %s"
                    % str(err).splitlines()[0][:200])


@pytest.fixture
def numpy_backend_restored(jax_device):
    yield
    bk.set_backend("numpy")


def test_saturated_task_kondo_spectrum_on_jax_matches_ed(numpy_backend_restored):
    bk.set_backend("jax")
    sc = impurity_chain(nb=3)
    _, ed = sc.get_kondo_spectrum(EV, mode="ED", **KW)
    _, tk = sc.get_kondo_spectrum(EV, mode="DMRG", submode="TASK", nkry=40, **KW)
    assert np.max(np.abs(tk - ed)) < 1e-9 * np.max(np.abs(ed))


def test_task_krylov_basis_stays_on_the_backend(numpy_backend_restored):
    """The basis, the environments and the matvec results are JAX arrays;
    only scalars and the Krylov-space matrices come back to the host."""
    bk.set_backend("jax")
    sc = impurity_chain(nb=3)
    sc.gs_energy()
    s = sc._session
    from dmrgpy.pyitensor.chain import _mps_arrays_lpr
    ts = task.TangentSpace(_mps_arrays_lpr(s.wf0), task.mpo_arrays(s.H))
    b = ts.project_operator(task.mpo_arrays(s._mpo(sc.Sz[0].to_terms())))
    store, n, Hq = task.block_krylov(ts, [b], 5)
    assert isinstance(store, jax.Array)
    assert isinstance(ts.LH[2], jax.Array)
    assert isinstance(Hq, np.ndarray) and Hq.shape == (n, n)


def _tangent_space(arrays):
    """The tangent space of the exported state and the three spin starts
    of site 0, on the backend that is set."""
    Ms, Ws, Wops = arrays
    ts = task.TangentSpace(Ms, Ws)
    return ts, [ts.project_operator(W) for W in Wops]


def test_the_residual_is_formed_on_the_backend(numpy_backend_restored,
                                               monkeypatch):
    """A Krylov space that is not exhausted (nkry=2: 9 of the 23 states)
    leaves a residual, whose Gram matrix is formed on the backend; only
    it, block by block, and vectors of the length of the basis come back,
    never a tangent vector."""
    sc = impurity_chain(nb=3)
    sc.gs_energy()
    s = sc._session
    from dmrgpy.pyitensor.chain import _mps_arrays_lpr
    arrays = (_mps_arrays_lpr(s.wf0), task.mpo_arrays(s.H),
              [task.mpo_arrays(s._mpo(getattr(sc, o)[0].to_terms()))
               for o in ("Sx", "Sy", "Sz")])
    ts, starts = _tangent_space(arrays)
    ref = task.block_krylov(ts, starts, 2, residual=True)
    bk.set_backend("jax")
    ts, starts = _tangent_space(arrays)
    N = ts.flat_size()
    shapes = []
    to_host = bk.to_host

    def spy(a):
        shapes.append(np.shape(a))
        return to_host(a)

    monkeypatch.setattr(bk, "to_host", spy)
    store, n, Hq, (j0, R) = task.block_krylov(ts, starts, 2, residual=True)
    monkeypatch.setattr(bk, "to_host", to_host)
    assert n == ref[1] < ts.dimension()
    assert R.shape == ref[3][1].shape and R.size and np.max(np.abs(R)) > 0
    assert not any(N in s for s in shapes)
    assert np.max(np.abs(R - ref[3][1])) < 1e-10 * np.max(np.abs(ref[3][1]))
    assert np.max(np.abs(Hq - ref[2])) < 1e-10 * np.max(np.abs(ref[2]))
