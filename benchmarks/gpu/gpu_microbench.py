"""Phase 0 of the pyitensor GPU port (docs/pyitensor_gpu_port_plan.md):
does a GPU beat a BLAS-pinned CPU on the *primitives* pyitensor actually
spends its time in, and above which bond dimension?

Standalone on purpose -- imports no dmrgpy, so it can be copied to a
cluster scratch directory on its own and cannot be confounded by anything
else in the package. Everything is complex128, because pyitensor's ITensor
is unconditionally complex128 (see pyitensor/tensor.py's docstring); that
means ZGEMM, i.e. ~4x the real-FLOP count of a same-shape real matmul, so
the crossover sits at a higher bond dimension than an FP32 ML workload
would suggest.

The four benchmarked primitives, and why exactly these:

* ``matvec`` -- the two-site effective-Hamiltonian application at the
  heart of DMRG's Lanczos solve and TDVP's Krylov step, i.e. the chain
  L * theta * W1 * W2 * R contracted in the same order (and therefore at
  the same cost) as pyitensor/dmrg.py's two_site_heff() and the plan that
  pyitensor/kernels.py precomputes for it. Dominated by two O(chi^3 d^2 w)
  tensordots. This is the op that decides the whole question: it is called
  dozens of times per bond, thousands of times per TDVP/METTS run.
* ``eigh`` -- the Gram route in pyitensor/svd.py's _svd_truncated(), which
  eigendecomposes the smaller of M M^dag / M^dag M rather than SVD-ing M.
* ``svd`` -- the exact-SVD fallback that same function takes when the
  smallest kept singular value is too small for the squared spectrum to
  resolve.
* ``qr`` -- used by svd.py's orthogonalization path.

plus ``transfer``, which is not an op at all but the thing that decides
whether a naive port helps: the H2D+D2H round trip for one theta-sized
array. If a matvec costs less than its own transfer, then any design that
converts per call (which is what pyitensor/kernels.py's existing JAX path
does today) cannot win no matter how fast the device is -- hence the
port's "arrays stay resident, only scalars come back" principle.

Shapes follow the DMRG two-site update: bond dimension chi, physical
dimension d (2 for S=1/2), MPO bond dimension w (5 for a Heisenberg
chain).
"""

import argparse
import json
import platform
import subprocess
import sys
import time


# ---------------------------------------------------------------- backends

class NumpyBackend:
    name = "numpy"

    def __init__(self):
        import numpy as np
        self.xp = np
        self.device = "cpu:%s" % platform.processor()

    def asarray(self, a):
        return self.xp.asarray(a, dtype=complex)

    def sync(self, a):
        """No-op: NumPy is synchronous. Returned so callers can use the
        same `sync(op())` shape for every backend."""
        return a

    def to_host(self, a):
        return self.xp.asarray(a)

    def compile(self, fn):
        """Hook for backends that can compile a kernel ahead of time; the
        default is to run the plain Python closure."""
        return fn

    def free_bytes(self):
        return None


class JaxBackend:
    name = "jax"

    def __init__(self):
        # x64 must be enabled before the first array is made, or every
        # array silently becomes complex64 and the comparison is
        # meaningless (pyitensor/kernels.py does the same thing).
        import jax
        jax.config.update("jax_enable_x64", True)
        import jax.numpy as jnp
        self.jax = jax
        self.xp = jnp
        self.device = ", ".join(str(d) for d in jax.devices())

    def asarray(self, a):
        return self.xp.asarray(a, dtype=complex)

    def sync(self, a):
        """JAX dispatch is asynchronous -- without this every timing here
        would measure enqueue time, not compute time."""
        return self.jax.block_until_ready(a)

    def to_host(self, a):
        import numpy as np
        return np.asarray(a)

    def compile(self, fn):
        return fn

    def free_bytes(self):
        try:
            stats = self.jax.devices()[0].memory_stats()
            return stats.get("bytes_limit", None)
        except Exception:
            return None


class TorchBackend:
    """Kept as a fallback candidate: same shapes, same timing protocol, so
    if JAX disappoints the comparison is already there. Its namespace is
    not NumPy-compatible enough for the `xp` shim, hence the aliases."""

    name = "torch"

    def __init__(self):
        import torch
        self.torch = torch
        self.dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.device = (torch.cuda.get_device_name(0)
                       if torch.cuda.is_available() else "cpu")

        class _XP:
            # torch spells tensordot's contracted axes `dims`, not `axes`
            # -- the third small namespace divergence found here, and a
            # concrete sample of what a torch `xp` adapter has to absorb.
            tensordot = staticmethod(
                lambda a, b, axes: torch.tensordot(a, b, dims=axes))

            class linalg:
                eigh = staticmethod(torch.linalg.eigh)
                svd = staticmethod(
                    lambda m: torch.linalg.svd(m, full_matrices=False))
                qr = staticmethod(torch.linalg.qr)

        self.xp = _XP()

    def asarray(self, a):
        return self.torch.as_tensor(
            a, dtype=self.torch.complex128, device=self.dev)

    def sync(self, a):
        if self.dev.type == "cuda":
            self.torch.cuda.synchronize()
        return a

    def to_host(self, a):
        # .cpu() first: numpy refuses a cuda tensor outright.
        return a.cpu().numpy()

    def compile(self, fn):
        return fn

    def free_bytes(self):
        if self.dev.type != "cuda":
            return None
        return self.torch.cuda.get_device_properties(0).total_memory


class JaxJitBackend(JaxBackend):
    """JAX with the matvec kernel jax.jit-compiled, which is the thing that
    decides whether eager jnp's per-call dispatch floor (measured at ~0.35 ms
    on an H200, vs torch's ~0.09 ms) is a property of JAX or only of running
    it eagerly. pyitensor/kernels.py already compiles this contraction as a
    single fused einsum for exactly this reason, so a jitted matvec is the
    realistic JAX design, not a trick.

    The kernel takes its arrays as *arguments* rather than closing over
    them: with them as closure constants XLA is free to constant-fold the
    whole contraction at compile time and the benchmark would measure
    nothing.
    """

    name = "jaxjit"

    def compile(self, fn):
        return self.jax.jit(fn)


BACKENDS = {"numpy": NumpyBackend, "jax": JaxBackend, "torch": TorchBackend,
            "jaxjit": JaxJitBackend}


# ------------------------------------------------------------------- ops

def _random_complex(rng, shape):
    import numpy as np
    return (rng.standard_normal(shape) + 1j * rng.standard_normal(shape)
            ).astype(np.complex128)


def make_matvec(be, chi, d, w, rng):
    """L * theta * W1 * W2 * R, same contraction order as
    pyitensor/dmrg.py's two_site_heff().

    Index names: a/b ket-and-bra left bonds, r right bond, m MPO bonds,
    s physical. Dominant cost is the first and last tensordot, both
    O(chi^3 d^2 w) -- exactly what a real two-site Lanczos iteration pays.
    """
    xp = be.xp
    L = be.asarray(_random_complex(rng, (chi, w, chi)))
    R = be.asarray(_random_complex(rng, (chi, w, chi)))
    W1 = be.asarray(_random_complex(rng, (w, d, d, w)))
    W2 = be.asarray(_random_complex(rng, (w, d, d, w)))
    theta = be.asarray(_random_complex(rng, (chi, d, d, chi)))

    def kernel(L, theta, W1, W2, R):
        # (a,m,b) x (a,s1,s2,r) -> (m,b,s1,s2,r)
        x = xp.tensordot(L, theta, axes=([0], [0]))
        # contract m and s1 with W1 (m,s1,s1',m') -> (b,s2,r,s1',m')
        x = xp.tensordot(x, W1, axes=([0, 2], [0, 1]))
        # contract m' and s2 with W2 (m',s2,s2',m'') -> (b,r,s1',s2',m'')
        x = xp.tensordot(x, W2, axes=([4, 1], [0, 1]))
        # contract r and m'' with R (r,m'',r') -> (b,s1',s2',r')
        return xp.tensordot(x, R, axes=([1, 4], [0, 1]))

    compiled = be.compile(kernel)

    def run():
        return compiled(L, theta, W1, W2, R)

    return run, {"theta_mib": (chi * d * d * chi * 16) / 2**20}


def make_eigh(be, chi, d, w, rng):
    """Hermitian eigendecomposition of a Gram matrix, the route
    pyitensor/svd.py's _svd_truncated() prefers over a full SVD."""
    xp = be.xp
    m = chi * d
    A = _random_complex(rng, (m, m))
    A = A + A.conj().T
    gram = be.asarray(A)

    def run():
        return xp.linalg.eigh(gram)

    return run, {"shape": [m, m]}


def make_svd(be, chi, d, w, rng):
    """Thin SVD of the two-site matrix, _svd_truncated()'s exact
    fallback."""
    xp = be.xp
    m, n = chi * d, d * chi
    mat = be.asarray(_random_complex(rng, (m, n)))

    def run():
        return xp.linalg.svd(mat)

    return run, {"shape": [m, n]}


def make_qr(be, chi, d, w, rng):
    xp = be.xp
    m, n = chi * d, d * chi
    mat = be.asarray(_random_complex(rng, (m, n)))

    def run():
        return xp.linalg.qr(mat)

    return run, {"shape": [m, n]}


def make_transfer(be, chi, d, w, rng):
    """One theta-sized host->device->host round trip. The number the
    port's whole design hinges on: compare it against `matvec`."""
    import numpy as np
    host = _random_complex(rng, (chi, d, d, chi))

    def run():
        dev = be.asarray(host)
        be.sync(dev)
        return be.to_host(dev)

    return run, {"mib": (chi * d * d * chi * 16) / 2**20}


OPS = {"matvec": make_matvec, "eigh": make_eigh, "svd": make_svd,
       "qr": make_qr, "transfer": make_transfer}


# ------------------------------------------------------------------ timing

def time_op(be, build, chi, d, w, rng, reps, warmup, budget):
    """Median of `reps` timed calls after `warmup` untimed ones.

    Every call is followed by the backend's own sync, so asynchronous
    dispatch (JAX, CUDA) is measured as compute rather than as enqueue.
    If a single call already exceeds `budget` seconds, that one call is
    the measurement -- so a CPU baseline at large chi degrades to one slow
    sample instead of blowing the job's walltime.
    """
    try:
        run, meta = build(be, chi, d, w, rng)
    except Exception as exc:                       # OOM, unsupported dtype
        return {"error": "%s: %s" % (type(exc).__name__, exc)}

    try:
        for _ in range(warmup):
            be.sync(run())

        t0 = time.perf_counter()
        be.sync(run())
        first = time.perf_counter() - t0

        if first > budget:
            times = [first]
        else:
            times = []
            for _ in range(reps):
                t0 = time.perf_counter()
                be.sync(run())
                times.append(time.perf_counter() - t0)
    except Exception as exc:
        return {"error": "%s: %s" % (type(exc).__name__, exc)}

    times.sort()
    out = {"median_s": times[len(times) // 2], "min_s": times[0],
           "n": len(times)}
    out.update(meta)
    return out


def gpu_info():
    try:
        out = subprocess.run(
            ["nvidia-smi", "--query-gpu=name,memory.total,driver_version",
             "--format=csv,noheader"],
            capture_output=True, text=True, timeout=30)
        return out.stdout.strip() or out.stderr.strip()
    except Exception as exc:
        return "nvidia-smi unavailable: %s" % exc


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--backends", default="numpy",
                    help="comma-separated: numpy, jax, torch")
    ap.add_argument("--chi", default="64,128,256,512,1024,2048")
    ap.add_argument("--d", type=int, default=2, help="physical dimension")
    ap.add_argument("--w", type=int, default=5, help="MPO bond dimension")
    ap.add_argument("--ops", default=",".join(OPS))
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--warmup", type=int, default=2)
    ap.add_argument("--budget", type=float, default=20.0,
                    help="seconds above which one call is the measurement")
    ap.add_argument("--json", default=None, help="write raw results here")
    args = ap.parse_args()

    chis = [int(c) for c in args.chi.split(",")]
    ops = [o for o in args.ops.split(",") if o]
    for o in ops:
        if o not in OPS:
            sys.exit("unknown op %r (have: %s)" % (o, ", ".join(OPS)))

    results = {"meta": {"python": sys.version.split()[0],
                        "platform": platform.platform(),
                        "nvidia_smi": gpu_info(),
                        "d": args.d, "w": args.w, "chi": chis,
                        "reps": args.reps, "warmup": args.warmup},
               "data": {}}

    for bname in args.backends.split(","):
        bname = bname.strip()
        if not bname:
            continue
        if bname not in BACKENDS:
            sys.exit("unknown backend %r" % bname)
        try:
            be = BACKENDS[bname]()
        except Exception as exc:
            print("[%s] unavailable: %s: %s" % (bname, type(exc).__name__, exc))
            results["data"][bname] = {"unavailable": str(exc)}
            continue

        print("\n=== %s === device: %s" % (bname, be.device))
        results["data"][bname] = {"device": be.device,
                                  "device_bytes": be.free_bytes(), "ops": {}}
        import numpy as np
        for op in ops:
            results["data"][bname]["ops"][op] = {}
            for chi in chis:
                rng = np.random.default_rng(1234)   # same data every backend
                r = time_op(be, OPS[op], chi, args.d, args.w, rng,
                            args.reps, args.warmup, args.budget)
                results["data"][bname]["ops"][op][str(chi)] = r
                if "error" in r:
                    print("  %-9s chi=%-5d  %s" % (op, chi, r["error"]))
                else:
                    print("  %-9s chi=%-5d  %10.4f ms  (n=%d)"
                          % (op, chi, 1e3 * r["median_s"], r["n"]))
                sys.stdout.flush()

    if args.json:
        with open(args.json, "w") as f:
            json.dump(results, f, indent=1)
        print("\nwrote %s" % args.json)


if __name__ == "__main__":
    main()
