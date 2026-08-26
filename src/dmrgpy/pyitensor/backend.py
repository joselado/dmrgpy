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

The dispatch floor, and the two knobs that fight it
--------------------------------------------------
Eager array operations carry a per-call dispatch floor -- measured at
~0.35 ms for eager jnp on the H200, versus ~0.07 ms for the same kernel
under `jax.jit`. That floor is why the device cannot win at small bond
dimension no matter how fast it is: a KPM run at kpmmaxm=40 issues ~20k
operations, so ~7 s of pure dispatch against a 5.8 s CPU baseline. It
amortizes away as chi grows (at kpmmaxm=160 the same dispatch is ~7% of a
103 s CPU run), which is the large-bond-dimension story the rest of the
port tells.

Two knobs lower that floor, and they only work together:

* `set_pad_bonds(K)` freezes every MPS bond at K, so the engine stops
  minting a fresh array shape every time a bond dimension changes;
* `set_jit()` + `jit()` fuse the hot composites -- a contraction's
  transpose+reshape+matmul, svd()'s Gram+eigh, a Lanczos reorthogonal-
  ization -- into a single XLA kernel each, so one dispatch replaces
  three to thirty.

`jax.jit` retraces per input shape (the same trap kernels.py's numba
section describes with the same cause), which is exactly what padding
removes; and padding on its own only reduces how many *eager* kernels get
compiled, it does not reduce how many get *dispatched*. Hence the default
`set_jit("auto")`: jit turns itself on precisely when `set_pad_bonds` has
made it safe, and stays off otherwise. Either can be forced.

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
    _JIT_CACHE.clear()   # jitted wrappers are bound to one array library
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


_JIT_MODE = "auto"
_JIT_CACHE = {}


def set_jit(mode="auto"):
    """Whether the hot composites registered through `jit()` below are
    compiled with `jax.jit` (True), left eager (False), or decided by
    `set_pad_bonds` ("auto", the default).

    Why "auto" is tied to padding. Jitting a composite replaces three to
    thirty eager dispatches with one, which is the single biggest lever on
    the ~0.35 ms/call floor that keeps a device from winning below
    chi ~ 120-160. But `jax.jit` traces once per distinct input *shape*,
    and DMRG mints a new shape every time a bond dimension changes -- so
    without `set_pad_bonds` the compile cost can exceed everything the
    fusion saves (measured on CPU at n=6/maxm=20: 672 compilations, 18.4 s
    of a 29.1 s run, already in eager mode). Padding collapses that shape
    zoo to one entry per operation, which is exactly the precondition jit
    needs; neither knob is worth much without the other.

    Forcing it is legitimate in both directions: True for a long run at
    fixed bond dimension that never padded (the shapes are stable anyway),
    False to isolate a compile-time effect while benchmarking.

    Never applies to the NumPy backend: `jit()` hands back the plain
    Python function there, so the host path is unchanged and costs one
    global comparison per call."""
    global _JIT_MODE
    if mode not in (True, False, "auto"):
        raise ValueError("set_jit: mode must be True, False or \"auto\", got %r" % (mode,))
    _JIT_MODE = mode
    return _JIT_MODE


def jit_mode():
    return _JIT_MODE


def jit_enabled():
    """Whether `jit()`-registered composites are currently compiled."""
    if _NAME != "jax":
        return False
    if _JIT_MODE == "auto":
        return _PAD_BONDS is not None
    return bool(_JIT_MODE)


def jit(fn, static_argnums=()):
    """Register `fn` as a hot composite: a wrapper that runs it under
    `jax.jit` when jitting is enabled, and calls it directly otherwise.

    The decision is made per call rather than at import, because the
    backend, the padding and the jit mode are all process-wide state a
    caller may change between two chains. That costs one string comparison
    per call -- nothing against a contraction, and *nothing at all* on
    NumPy, which is the path that must not regress.

    `static_argnums` marks the arguments that are shapes/flags rather than
    arrays, so XLA specializes on them instead of trying to trace them.
    They must be hashable (tuples, not lists)."""
    key = (fn, static_argnums)

    def wrapper(*args):
        if _NAME == "jax" and jit_enabled():
            compiled = _JIT_CACHE.get(key)
            if compiled is None:
                compiled = _JAX.jit(fn, static_argnums=static_argnums)
                _JIT_CACHE[key] = compiled
            return compiled(*args)
        return fn(*args)

    wrapper.__name__ = getattr(fn, "__name__", "jitted")
    wrapper.__doc__ = fn.__doc__
    return wrapper


def compilations():
    """How many distinct (function, shape) kernels the jitted composites
    have traced so far -- the number `set_pad_bonds` exists to collapse.
    Returns None when jitting is off. Diagnostic only; benchmarks/gpu/
    reads it to report compile pressure alongside wall time."""
    if not jit_enabled():
        return None
    total = 0
    for compiled in _JIT_CACHE.values():
        try:
            total += compiled._cache_size()
        except Exception:      # private API, version-dependent
            return None
    return total


def sync(a):
    """Block until `a` is actually computed. JAX dispatch is asynchronous,
    so any timing that does not call this measures enqueue, not compute."""
    if _NAME == "jax" and a is not None:
        _JAX.block_until_ready(a)
    return a
