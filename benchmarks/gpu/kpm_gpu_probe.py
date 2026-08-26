"""Is the KPM dynamical correlator worth porting to GPU, and which part of
it? Companion to gpu_microbench.py/pyitensor_gpu_probe.py -- see
docs/pyitensor_gpu_port_plan.md.

The KPM correlator has a *different* cost profile from ground-state DMRG,
which is why it needs its own measurement rather than inheriting the
ground-state conclusion. DMRG's sweep is dominated by the two-site
effective-Hamiltonian matvec (a big dense contraction, the primitive the
GPU wins most decisively). The KPM moment recursion in
pyitensor/chain.py's _kpm_moments_full/_kpm_moments_accelerated is instead
a long serial chain of

    a' = 2 * (M a) - a_prev,      truncated to kpmmaxm after every step

so each of its ~n moments costs one MPO application plus one MPS sum, and
*both* of those truncate through svd.py -- the primitive the GPU is
*worst* at (measured slower than CPU below chi=256 on an H200). Profiled
on CPU at n=12, kpmmaxm=40, ~150 moments: 78% of the run is in the moment
recursion, 49% inside svd(), 25% inside the Gram-route eigh alone, and
19% in tensordot spread over 17711 calls. So the KPM port is a
factorization-and-many-small-calls problem, not a big-contraction problem.

What this script measures, per (n, kpmmaxm):

* wall time split between the ground state (a DMRG sweep, matvec-bound)
  and the moment recursion (factorization-bound), since only the second
  is what "porting KPM" means;
* the recursion's own internal split -- svd / eigh / tensordot / everything
  else -- by profiling the call rather than guessing, so the port can
  target the dominant term;
* how all of that scales with kpmmaxm, which is the axis that decides the
  question: at typical kpmmaxm (30-60) the matrices are far below the
  crossover the primitive benchmark found, while at kpmmaxm >= 256 the
  Gram eigh is 4-27x faster on the device.

Read the BUCKETS comment below before trusting any percentage this
prints: cProfile can only attribute *callables*, and it inflates this
particular workload by ~1.7x because the workload is millions of small
calls. Both effects push the apparent Python-side share up, which is how
an earlier version of this script concluded (wrongly) that KPM was
half host-side overhead.

Correctness is checked at every point via the zeroth-moment sum rule
(int S_AA domega == <GS|A A|GS> == 1/4 for A=Sz on spin-1/2), the same
invariant tests/test_kpm_dynamical_correlator_python.py asserts, so a
timing point can never be a silently wrong spectrum.
"""

import argparse
import cProfile
import json
import os
import pstats
import sys
import time


def build_chain(n, maxm, kpmmaxm, nsweeps):
    from dmrgpy import spinchain
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


# Profile buckets, and an important caveat about what they can and cannot
# show. cProfile only sees *callables*: named functions like np.linalg.eigh
# and np.tensordot appear as their own entries, but arithmetic invoked
# through operators or methods -- `am @ bm`, `mat.conj().T`, `np.sqrt`,
# fancy slicing -- is executed in C without a Python frame and is charged
# to the *self time of the calling function*. So svd.py's or tensor.py's
# self time is a mixture of real BLAS work and Python bookkeeping, and is
# NOT attributable to either.
#
# This bit an earlier version of this script, which split the run into
# "array math" and "host bookkeeping" and concluded that 55-63% of KPM was
# host-side overhead with an Amdahl ceiling of ~1.8x on any GPU port. That
# was wrong: most of svd.py's self time is the Gram matmul and the
# recovery matmul, i.e. exactly the arithmetic a device would accelerate.
# The A/B in this repo's history settles it from the other direction --
# replacing np.tensordot with an inline transpose+reshape+matmul in
# ITensor.__mul__, which removes tensordot's Python-level overhead
# entirely, moved the end-to-end KPM time by only a few percent, which it
# could not have done if Python overhead were half the run.
#
# So the buckets below are reported as three honest groups instead:
#   named_linalg     -- a hard LOWER bound on array arithmetic
#   mixed_self       -- BLAS + bookkeeping, not separable by cProfile
#   bookkeeping      -- a hard lower bound on Python-only overhead
# and no Amdahl ceiling is claimed. The array-math share is somewhere
# between named_linalg and named_linalg+mixed_self; pinning it down needs
# an end-to-end device measurement (Phase 2), not a profile.
BUCKETS = [
    # --- hard lower bound on array arithmetic: named numpy entry points ---
    ("eigh", lambda f, fn: fn == "eigh"),
    ("svd_qr", lambda f, fn: fn in ("svd", "qr", "eigvalsh", "eigh_lapack")),
    ("tensordot", lambda f, fn: fn == "tensordot"),
    ("numpy_other", lambda f, fn: "numpy" in f or "scipy" in f),
    # --- inseparable mixture of inline BLAS and Python bookkeeping ---
    ("svd_py", lambda f, fn: f.endswith("svd.py")),
    ("tensor_py", lambda f, fn: f.endswith("tensor.py")),
    ("mps_algebra", lambda f, fn: f.endswith(("mpsalgebra.py", "mpscontainer.py"))),
    ("chain_kpm", lambda f, fn: f.endswith(("chain.py", "dmrg.py", "kernels.py"))),
    # --- hard lower bound on Python-only overhead ---
    ("bookkeeping", lambda f, fn: f.endswith("index.py")),
]

NAMED_LINALG = ("eigh", "svd_qr", "tensordot", "numpy_other")
MIXED = ("svd_py", "tensor_py", "mps_algebra", "chain_kpm")


def bucket_profile(prof):
    """Self time per bucket. Self time (not cumulative) so the buckets
    partition the run instead of double-counting a call tree."""
    stats = pstats.Stats(prof)
    out = {name: 0.0 for name, _ in BUCKETS}
    out["other"] = 0.0
    for (fname, _lineno, funcname), (_cc, _nc, tottime, _ct, _cs) in stats.stats.items():
        label = "other"
        for name, match in BUCKETS:
            if match(fname, funcname):
                label = name
                break
        out[label] += tottime
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--src", default=None,
                    help="path to dmrgpy's src/ (prepended to sys.path)")
    ap.add_argument("--n", type=int, default=12, help="chain length")
    ap.add_argument("--kpmmaxm", default="20,40,80,160",
                    help="comma-separated KPM bond dimensions to sweep")
    ap.add_argument("--maxm", type=int, default=30,
                    help="ground-state bond dimension")
    ap.add_argument("--nsweeps", type=int, default=6)
    ap.add_argument("--delta", type=float, default=0.08,
                    help="KPM resolution; sets the moment count")
    ap.add_argument("--json", default=None)
    ap.add_argument("--plot", default=None, help="write a PNG here")
    args = ap.parse_args()

    if args.src:
        sys.path.insert(0, os.path.abspath(args.src))
    import numpy as np

    out = {"meta": {"n": args.n, "maxm": args.maxm, "delta": args.delta,
                    "nsweeps": args.nsweeps, "python": sys.version.split()[0]},
           "runs": []}
    try:
        import jax
        out["meta"]["jax_devices"] = [str(d) for d in jax.devices()]
    except Exception as exc:
        out["meta"]["jax_error"] = str(exc)
    from dmrgpy.pyitensor import kernels
    out["meta"]["use_jax"] = bool(kernels.USE_JAX)

    es = np.linspace(-1.0, 6.0, 400)
    es_wide = np.linspace(-8.0, 8.0, 2000)

    for kpmmaxm in [int(k) for k in args.kpmmaxm.split(",")]:
        sc = build_chain(args.n, args.maxm, kpmmaxm, args.nsweeps)

        t0 = time.perf_counter()
        sc.gs_energy(mode="DMRG")
        t_gs = time.perf_counter() - t0

        # Two runs of the same call: one clean (the only trustworthy wall
        # time) and one profiled (for the cost split). cProfile charges
        # per-Python-call overhead, and this workload makes millions of
        # small calls, so the profiled run is substantially slower AND
        # skewed toward call-heavy buckets -- measured 1.71x inflation at
        # kpmmaxm=160 (102.7s clean vs 175.9s profiled). Reporting only the
        # profiled number, as an earlier version did, both overstates the
        # runtime and overstates the Python-side share of it.
        t0 = time.perf_counter()
        x, y = sc.get_dynamical_correlator(mode="DMRG", submode="KPM",
                                            name=(sc.Sz[0], sc.Sz[0]),
                                            es=es, delta=args.delta)
        t_kpm = time.perf_counter() - t0

        prof = cProfile.Profile()
        prof.enable()
        t0 = time.perf_counter()
        sc.get_dynamical_correlator(mode="DMRG", submode="KPM",
                                     name=(sc.Sz[0], sc.Sz[0]),
                                     es=es, delta=args.delta)
        t_profiled = time.perf_counter() - t0
        prof.disable()

        # Sum rule on a wide window: no timing point is trusted without it.
        xw, yw = sc.get_dynamical_correlator(mode="DMRG", submode="KPM",
                                              name=(sc.Sz[0], sc.Sz[0]),
                                              es=es_wide, delta=args.delta)
        weight = float(np.trapezoid(np.real(yw), xw))

        row = {"kpmmaxm": kpmmaxm, "t_groundstate": t_gs, "t_kpm": t_kpm,
               "t_profiled": t_profiled,
               "profiler_inflation": t_profiled / t_kpm,
               "buckets": bucket_profile(prof), "sum_rule": weight,
               "sum_rule_error": abs(weight - 0.25),
               "peak": float(xw[np.argmax(np.real(yw))])}
        out["runs"].append(row)

        b = row["buckets"]
        named = sum(b[k] for k in NAMED_LINALG)
        mixed = sum(b[k] for k in MIXED)
        # of the profiled total, since that is what the buckets partition
        row["named_linalg_fraction"] = named / t_profiled
        row["mixed_fraction"] = mixed / t_profiled
        row["bookkeeping_fraction"] = b["bookkeeping"] / t_profiled
        print("kpmmaxm=%-5d gs %6.2fs  kpm %7.2fs (profiled %7.2fs, "
              "%.2fx inflation)   of profiled: named linalg %3.0f%% "
              "(eigh %2.0f svd/qr %2.0f tdot %2.0f)  mixed %3.0f%%  "
              "bookkeeping %3.0f%%   sum rule %.6f (err %.1e)"
              % (kpmmaxm, t_gs, t_kpm, t_profiled, t_profiled / t_kpm,
                 100 * named / t_profiled,
                 100 * b["eigh"] / t_profiled, 100 * b["svd_qr"] / t_profiled,
                 100 * b["tensordot"] / t_profiled, 100 * mixed / t_profiled,
                 100 * b["bookkeeping"] / t_profiled,
                 weight, row["sum_rule_error"]))
        sys.stdout.flush()

    if args.json:
        with open(args.json, "w") as f:
            json.dump(out, f, indent=1)
        print("wrote %s" % args.json)

    if args.plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib unavailable; skipped the figure")
            return
        ks = [r["kpmmaxm"] for r in out["runs"]]
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
        axes[0].plot(ks, [r["t_kpm"] for r in out["runs"]], "o-",
                     label="KPM moment recursion")
        axes[0].plot(ks, [r["t_groundstate"] for r in out["runs"]], "s--",
                     label="ground state")
        axes[0].set_xscale("log", base=2)
        axes[0].set_yscale("log")
        axes[0].set_xlabel("kpmmaxm")
        axes[0].set_ylabel("wall time [s]")
        axes[0].set_title("where the time goes (n=%d)" % args.n)
        axes[0].legend()
        axes[0].grid(alpha=.3)

        labels = list(NAMED_LINALG) + list(MIXED) + ["bookkeeping", "other"]
        bottom = [0.0] * len(ks)
        for lab in labels:
            vals = [100 * r["buckets"][lab] / r["t_kpm"] for r in out["runs"]]
            axes[1].bar([str(k) for k in ks], vals, bottom=bottom, label=lab)
            bottom = [b + v for b, v in zip(bottom, vals)]
        axes[1].set_xlabel("kpmmaxm")
        axes[1].set_ylabel("% of KPM wall time (self time)")
        axes[1].set_title("cost split (see BUCKETS: the middle group is\n"
                          "inline BLAS + bookkeeping, not separable)")
        axes[1].legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(args.plot, dpi=140)
        print("wrote %s" % args.plot)


if __name__ == "__main__":
    main()
