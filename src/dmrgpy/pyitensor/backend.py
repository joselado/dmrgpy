"""Which array library pyitensor's ITensors are made of: NumPy (the
default, and the only one before this module existed) or JAX, so the whole
engine can run on a GPU.

Design, and why it is this small
-------------------------------
The port needed exactly one property to be cheap: **arrays must stay on
the device**. Phase 0 of docs/pyitensor_gpu_port_plan.md measured, on an
NVIDIA H200, that a single host->device->host round trip for one
theta-sized array costs *more* than the two-site matvec it would feed
(22.0 ms vs 8.3 ms at chi=1024). Anything that converts per call loses
however fast the device is -- which is precisely why the pre-existing
per-call JAX path in kernels.py made GPU runs 5-11x *slower*.

So there is exactly one conversion point in the whole engine:
ITensor.__init__ (via asarray below). Everything built from a converted
array stays on the device, because every operation pyitensor performs on
tensor data is either an ndarray *method* (`.transpose`, `.reshape`,
`.conj`, `@`, arithmetic operators) -- which NumPy and JAX both provide
with identical semantics -- or one of the handful of free functions
routed through this module. That is why the diff for this port is small:
most of tensor.py never needed a namespace at all.

Host round trips are confined to three places, all of them O(1) or
O(chi) rather than O(chi^2):

* `scalar()` -- an energy, an overlap, a Lanczos alpha/beta;
* `to_host()` on a singular-value vector, so svd.py's truncation rule
  (a cumulative sum plus data-dependent branching, worthless on a device)
  runs on the host;
* measurement results on their way back to Python.

What this deliberately does not do
----------------------------------
No `jax.jit`. Every operation here runs eagerly, which on GPU carries a
per-call dispatch floor -- measured at ~0.35 ms for eager jnp on the H200,
versus ~0.07 ms for the same kernel under jit. That floor is why the
device cannot win at small bond dimension no matter what: a KPM run at
kpmmaxm=40 issues ~20k operations, so ~7 s of pure dispatch against a
5.8 s CPU baseline. It amortizes away as chi grows (at kpmmaxm=160 the
same dispatch is ~7% of a 103 s CPU run), which is the same
large-bond-dimension story the rest of the plan tells. Jitting the hot
composites is a possible later step; it is not needed to answer whether
the port pays off, and it interacts badly with DMRG's ever-changing bond
shapes (jax.jit retraces per shape -- see kernels.py's numba section for
the same trap with the same cause).

Usage
-----
    from dmrgpy.pyitensor import backend
    backend.set_backend("jax")      # or "numpy" (default)
    print(backend.device_info())

Switch before building a chain: existing ITensors keep whatever array
type they were made with, and mixing the two in one contraction is not
supported (it would silently transfer per call, the exact thing this
module exists to avoid).
"""

import numpy as np

_NAME = "numpy"
_XP = np
_JAX = None


def set_backend(name):
    """Select the array library for arrays created from now on. "numpy"
    (default) or "jax". Raises rather than falling back silently -- a
    GPU run that quietly became a CPU run is a benchmark that lies."""
    global _NAME, _XP, _JAX
    if name == "numpy":
        _NAME, _XP, _JAX = "numpy", np, None
        return _NAME
    if name == "jax":
        import jax
        # Must precede the first array: JAX silently makes everything
        # complex64 otherwise, and this engine is complex128 throughout.
        jax.config.update("jax_enable_x64", True)
        import jax.numpy as jnp
        _NAME, _XP, _JAX = "jax", jnp, jax
        return _NAME
    raise ValueError("unknown pyitensor backend %r (have: numpy, jax)" % name)


def backend_name():
    return _NAME


def xp():
    """The active array namespace. Resolve once per function, never per
    element -- this is a global lookup on a hot path otherwise."""
    return _XP


def is_device():
    """True when arrays live somewhere a host pointer cannot reach."""
    return _NAME != "numpy"


def device_info():
    if _NAME == "numpy":
        return "numpy (host)"
    return "jax: %s" % ", ".join(str(d) for d in _JAX.devices())


def asarray(a):
    """The single conversion point of the whole engine (ITensor.__init__).
    A host array becomes a device array here and stays there; a device
    array passes through untouched."""
    return _XP.asarray(a, dtype=complex)


def zeros(shape):
    return _XP.zeros(shape, dtype=complex)


def eye(n):
    return _XP.eye(n, dtype=complex)


def to_host(a):
    """Device -> host. Every call is a synchronization point, so keep them
    O(chi) or smaller; see this module's docstring."""
    return np.asarray(a)


def scalar(a):
    """One number back to Python (energy, overlap, Lanczos coefficient)."""
    return complex(np.asarray(a))


def setblock(arr, idx, value):
    """arr[idx] = value, functionally.

    NumPy assigns in place and returns the same buffer; JAX arrays are
    immutable, so this returns a *new* array via .at[].set(). Callers must
    therefore use the return value rather than relying on mutation --
    mpsalgebra.py's direct-sum construction is the only place in the
    engine that needs it."""
    if _NAME == "numpy":
        arr[idx] = value
        return arr
    return arr.at[idx].set(value)


_PAD_BONDS = None


def set_pad_bonds(dim):
    """Pad every MPS bond produced by svd() up to exactly `dim` with zero
    singular values (None = off, the default and the historical behavior).

    This exists for one reason: XLA compiles a kernel per distinct
    (operation, shape), and DMRG/KPM mint a new shape every time a bond
    dimension changes. Measured on CPU at n=6/maxm=20, a single
    ground-state solve triggered 672 compilations costing 18.4 s of a
    29.1 s run. Freezing every bond at one dimension collapses that shape
    zoo to a single entry per operation, so each kernel is compiled once
    and every later call reuses it.

    It is exact, not an approximation. The padding appends *zero*
    singular values and zero vectors after the truncation has already
    chosen what to keep, so the represented state is unchanged: the extra
    columns of U and rows of V contract to zero everywhere. It is not
    free, though -- near the chain edges the true bond dimension is
    2, 4, 8, ... and padding those up to `dim` does arithmetic on blocks
    that are known to be zero. On a device that trade is usually worth it
    (uniform shapes, one compile, better occupancy); on the host it is
    pure waste, which is why it is opt-in rather than automatic.

    Padding happens after the factorization, never before, so the Gram
    route's `Vh = (W^dag M) / S` division is only ever taken over the
    singular values that survived truncation -- forcing mindim up to
    `dim` instead would divide by ~0 and push every call into the exact-
    SVD fallback, which is the one primitive a GPU is worst at.
    """
    global _PAD_BONDS
    _PAD_BONDS = int(dim) if dim else None
    return _PAD_BONDS


def pad_bonds():
    return _PAD_BONDS


def sync(a):
    """Block until `a` is actually computed. JAX dispatch is asynchronous,
    so any timing that does not call this measures enqueue, not compute."""
    if _NAME == "jax" and a is not None:
        _JAX.block_until_ready(a)
    return a
