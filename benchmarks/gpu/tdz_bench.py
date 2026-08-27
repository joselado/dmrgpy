"""The complex-time-evolution dynamical correlator (submode="TDZ",
arXiv:2311.10909) on host versus device, with each optimization
attributable on its own.

What this measures, and why it is shaped this way
-------------------------------------------------

TDZ is one long chain of *small* array operations: a two-site TDVP step
per time step (2 half-sweeps x N bonds x a Krylov exponentiation per
bond), plus n_max+1 full-chain overlaps per step to get the phi^(n) the
real-axis reconstruction needs. Nothing in it is a big GEMM, so on a
device its cost is made almost entirely of per-call dispatch and of
host synchronizations that stall JAX's asynchronous queue. Three changes
attack exactly that:

  krylov  tdvp.set_krylov_defer_sync -- the Krylov exponentiator keeps
          alpha/beta on the device and brings a whole block home at once,
          instead of two `float(to_host(...))` synchronizations per
          Lanczos iteration (~1000 per time step at n=30). It speculates
          past its own stopping point and rolls back to the exact same
          Krylov dimension, so the answer is unchanged.
  bras    tdz.set_batched_bras -- the n_max+1 phi^(n) overlaps as one
          batched sweep (mpsalgebra.BatchedBras) instead of n_max+1
          separate ones: same arithmetic, a fifth of the dispatches.
  (always on) svd()'s diagonal S tensor is now built from the
          device-resident spectrum rather than np.diag() of its host
          copy, which removed an O(chi^2) host->device transfer per
          factorization -- i.e. per bond, per half-sweep, per step. That
          one is not a trade-off and has no switch; it is included in
          every row below.

The four (krylov, bras) combinations are run in one process so a single
job yields the whole attribution table, and both are *scheduling*
changes, so every configuration must return the same spectrum -- which
is checked here rather than assumed, and reported as a max deviation
against the baseline row.

Cold versus warm is reported for the same reason port_speedup.py reports
it: XLA compiles a kernel per distinct (operation, shape), so a
one-shot script and a long-running process see very different costs.

The CPU baseline must run as its own job on a CPU partition (CPU-only
workloads are not permitted on the cluster's GPU partitions), so the two
halves are combined afterwards with --from-json, exactly as the rest of
this directory does.

Expectation worth stating up front, so the result is read honestly: TDZ's
whole selling point is that a complex-time contour damps high-energy
content and so needs a *small* bond dimension (the paper reports chi ~
20-30 against 500-700 for real-time evolution), while this port's
measured crossover is chi ~ 120-160. At its natural operating point TDZ
is therefore on the wrong side of that line, and the maxm sweep here is
the deliverable rather than a headline speedup: it says where the
crossover for *this* calculation actually falls.
"""

import argparse
import json
import os
import sys
import time


CONFIGS = [
    ("base", False, False),
    ("krylov", True, False),
    ("bras", False, True),
    ("both", True, True),
]


def _exchange(sc, i, j):
    return sc.Sx[i] * sc.Sx[j] + sc.Sy[i] * sc.Sy[j] + sc.Sz[i] * sc.Sz[j]


def build_chain(n, maxm, nsweeps, seed):
    """A uniform S=1/2 Heisenberg chain -- the same model the rest of this
    directory's dynamical-correlator benchmarks use, and the right one
    here: unlike the ground-state sweeps (see port_speedup.build_chain on
    why a 1D chain cannot show a ground-state GPU win), a *time-evolved*
    state grows entanglement with time, so the bond dimension actually
    reaches the requested maxm rather than saturating at the ground
    state's own chi ~ 60."""
    import numpy as np
    from dmrgpy import spinchain

    np.random.seed(seed)
    sc = spinchain.Spin_Chain(["S=1/2" for _ in range(n)],
                              itensor_version="python")
    h = 0
    for i in range(n - 1):
        h = h + _exchange(sc, i, i + 1)
    sc.set_hamiltonian(h)
    sc.maxm = maxm
    sc.nsweeps = nsweeps
    return sc


def time_tdz(n, maxm, nt, dt, nsweeps, seed, reps, alpha0, n_max):
    """(times, spectrum) for `reps` runs of one TDZ correlator.

    The ground state is solved inside the timed region on purpose: it is
    part of what a user pays for one correlator, and it is the same cost
    in every configuration, so it does not distort the comparison between
    them. It is reported separately as well, so a reader can subtract it.
    """
    import numpy as np
    from dmrgpy.pyitensor import backend as bk

    es = np.linspace(-0.5, 4.0, 200)
    times = []
    spectrum = None
    gs_time = None
    for rep in range(reps):
        sc = build_chain(n, maxm, nsweeps, seed)
        t0 = time.perf_counter()
        sc.get_gs()
        bk.sync(getattr(sc.wf0.cpp_handle.A(1), "array", None))
        t_gs = time.perf_counter()
        _es, spectrum = sc.get_dynamical_correlator(
                mode="DMRG", submode="TDZ", name=(sc.Sz[n // 2], sc.Sz[n // 2]),
                alpha0=alpha0, n_max=n_max, dt=dt, nt=nt, es=es)
        times.append(time.perf_counter() - t0)
        if gs_time is None:
            gs_time = t_gs - t0
    return times, np.asarray(spectrum), gs_time


def run_backend(name, args, results):
    import numpy as np
    from dmrgpy import tdz as tdzmod
    from dmrgpy.pyitensor import backend as bk
    from dmrgpy.pyitensor import tdvp

    try:
        bk.set_backend(name)
    except Exception as exc:
        print("[%s] unavailable: %s" % (name, exc))
        return
    print("\n=== backend %s: %s ===" % (name, bk.device_info()))
    results["meta"]["device_" + name] = bk.device_info()
    results["tdz"][name] = {}

    for maxm in [int(m) for m in args.maxm.split(",")]:
        if args.pad_bonds:
            bk.set_pad_bonds(maxm)
        results["tdz"][name][str(maxm)] = {}
        reference = None
        for label, defer, batched in CONFIGS:
            if label not in args.configs.split(","):
                continue
            # On NumPy the deferred-sync path is pure waste (to_host is a
            # no-op there), so it is only meaningful to force it on a
            # device; it is still run on the host so the table has the
            # same rows everywhere and any host regression would show.
            tdvp.set_krylov_defer_sync(defer)
            tdzmod.set_batched_bras(batched)
            try:
                times, spectrum, gs_time = time_tdz(
                        args.n, maxm, args.nt, args.dt, args.nsweeps,
                        args.seed, args.reps, args.alpha0, args.n_max)
            except Exception as exc:
                print("  maxm=%-4d %-7s FAILED %s: %s"
                      % (maxm, label, type(exc).__name__, exc))
                results["tdz"][name][str(maxm)][label] = {"error": str(exc)}
                continue
            warm = min(times[1:]) if len(times) > 1 else times[0]
            if reference is None:
                reference = spectrum
                dev = 0.0
            else:
                dev = float(np.max(np.abs(spectrum - reference)))
            entry = {"times": times, "cold": times[0], "warm": warm,
                     "gs_time": gs_time, "max_dev_vs_base": dev,
                     "compilations": bk.compilations()}
            results["tdz"][name][str(maxm)][label] = entry
            print("  maxm=%-4d %-7s cold %8.2fs  warm %8.2fs  "
                  "(gs %6.2fs)  dev %.2e  kernels %s"
                  % (maxm, label, times[0], warm, gs_time, dev,
                     entry["compilations"]))
            sys.stdout.flush()
        tdvp.set_krylov_defer_sync(None)
        tdzmod.set_batched_bras(True)
        bk.set_pad_bonds(None)


def report(results):
    """Host/device table per bond dimension, plus the attribution."""
    print("\n" + "=" * 78)
    print("TDZ dynamical correlator: n=%s, nt=%s, dt=%s, alpha0=%s, n_max=%s"
          % (results["meta"]["n"], results["meta"]["nt"], results["meta"]["dt"],
             results["meta"]["alpha0"], results["meta"]["n_max"]))
    print("=" * 78)
    tdz = results.get("tdz", {})
    maxms = sorted({m for b in tdz.values() for m in b}, key=int)
    for maxm in maxms:
        print("\nmaxm = %s" % maxm)
        print("  %-8s %-10s %12s %12s %10s" % ("backend", "config", "cold(s)",
                                               "warm(s)", "dev"))
        for name in sorted(tdz):
            for label, _d, _b in CONFIGS:
                e = tdz[name].get(maxm, {}).get(label)
                if not e:
                    continue
                if "error" in e:
                    print("  %-8s %-10s %s" % (name, label, e["error"][:50]))
                    continue
                print("  %-8s %-10s %12.2f %12.2f %10.1e"
                      % (name, label, e["cold"], e["warm"],
                         e["max_dev_vs_base"]))
        host = tdz.get("numpy", {}).get(maxm, {}).get("both")
        dev = tdz.get("jax", {}).get(maxm, {}).get("both")
        if host and dev and "error" not in host and "error" not in dev:
            print("  -> device speedup at this maxm: %.2fx cold, %.2fx warm"
                  % (host["cold"] / dev["cold"], host["warm"] / dev["warm"]))
        base = tdz.get("jax", {}).get(maxm, {}).get("base")
        if base and dev and "error" not in base and "error" not in dev:
            print("  -> optimizations on device: %.2fx cold, %.2fx warm"
                  % (base["cold"] / dev["cold"], base["warm"] / dev["warm"]))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--src", default=None,
                    help="checkout's src/ to put at the front of sys.path")
    ap.add_argument("--backends", default="numpy",
                    help="comma-separated: numpy, jax")
    ap.add_argument("--configs", default="base,krylov,bras,both",
                    help="which optimization combinations to time")
    ap.add_argument("--n", type=int, default=30, help="chain length")
    ap.add_argument("--maxm", default="30,60,120,240",
                    help="MPS bond dimensions to sweep")
    ap.add_argument("--nt", type=int, default=100,
                    help="complex-time steps (the cost is linear in this; "
                         "a production TDZ run uses ~1000, which is why the "
                         "per-step cost is what this reports)")
    ap.add_argument("--dt", type=float, default=0.1)
    ap.add_argument("--alpha0", type=float, default=0.1)
    ap.add_argument("--n-max", dest="n_max", type=int, default=4)
    ap.add_argument("--nsweeps", type=int, default=4)
    ap.add_argument("--reps", type=int, default=2,
                    help=">=2 so the cold (compiling) run and the warm runs "
                         "are separated")
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument("--pad-bonds", action="store_true",
                    help="freeze every MPS bond at the sweep's bond dimension "
                         "(backend.set_pad_bonds), which is also what turns "
                         "backend.set_jit(\"auto\") on")
    ap.add_argument("--json", default=None)
    ap.add_argument("--from-json", default=None, nargs="*",
                    help="merge previous runs' JSON and re-print the report "
                         "without recomputing (how the CPU and GPU jobs are "
                         "combined)")
    args = ap.parse_args()

    if args.src:
        sys.path.insert(0, os.path.abspath(args.src))

    if args.from_json:
        merged = None
        for path in args.from_json:
            with open(path) as fh:
                part = json.load(fh)
            if merged is None:
                merged = part
            else:
                for name, block in part.get("tdz", {}).items():
                    merged.setdefault("tdz", {}).setdefault(name, {}).update(block)
                merged["meta"].update(part.get("meta", {}))
        report(merged)
        return

    results = {"meta": {"n": args.n, "nt": args.nt, "dt": args.dt,
                        "alpha0": args.alpha0, "n_max": args.n_max,
                        "maxm": args.maxm, "nsweeps": args.nsweeps,
                        "reps": args.reps, "seed": args.seed,
                        "pad_bonds": bool(args.pad_bonds),
                        "python": sys.version.split()[0]},
               "tdz": {}}
    try:
        import subprocess
        results["meta"]["nvidia_smi"] = subprocess.run(
            ["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"],
            capture_output=True, text=True, timeout=30).stdout.strip()
    except Exception:
        pass

    for name in args.backends.split(","):
        name = name.strip()
        if name:
            run_backend(name, args, results)

    report(results)
    if args.json:
        os.makedirs(os.path.dirname(os.path.abspath(args.json)), exist_ok=True)
        with open(args.json, "w") as fh:
            json.dump(results, fh, indent=1)
        print("\nwrote %s" % args.json)


if __name__ == "__main__":
    main()
