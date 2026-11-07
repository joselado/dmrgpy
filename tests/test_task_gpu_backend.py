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
