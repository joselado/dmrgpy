"""METTS (arXiv:1002.1305) and dynamical METTS (arXiv:2405.18484) on host
versus device: where the crossover falls for a sampling method.

What this measures, and why it is shaped this way
-------------------------------------------------

METTS was the last big item on docs/pyitensor_gpu_port_plan.md Sec. 9's
list, and it is the one whose cost is almost entirely *inherited*: a
sample is one imaginary-time TDVP evolution (ported, item 2), the
observables are mpsalgebra.inner (ported), the truncations are svd.py
(ported). The only piece that was still host-bound is the collapse step,
and it is O(n) small operations per sample against thousands inside the
evolution. So the honest question here is not "how much did porting
collapse_to_cps() buy" -- it is "at what bond dimension is METTS worth a
device at all", which is a property of the calculation rather than of the
port.

Two calculations, because they sit on opposite sides of that question:

  vev  metts_vev: nwarmup+nsamples imaginary-time evolutions of a product
       state, each measured against an observable. METTS states are
       *minimally entangled* by construction -- that is the method's whole
       selling point -- so chi stays small and this is expected to be on
       the wrong side of the port's chi ~ 120-160 crossover. Worth
       measuring rather than assuming, exactly as TDZ was.
  dc   metts_dynamical_correlator: the same Markov chain, plus two
       *real-time* evolutions per retained sample (|v(t)>, |w(t)>) whose
       entanglement grows with t by construction, so chi climbs into the
       hundreds the way TDVP's does. This is where a device should win.

Both defaults below were set by measurement, not by analogy with the other
scripts here, and the first measurement was a negative one. At T=0.5 with
nt=20 (t <= 2) the whole sweep is *flat* in maxm -- 1.02 s for vev and
5.2 s for dc at both maxm=60 and maxm=240, on one core -- because neither
calculation's bond dimension ever reaches the cap, so a sweep run that way
measures nothing about bond dimension at all. What does bite, measured the
same way:

  dc, nt=100 (t <= 10):  21.6 s / 71.7 s / 1107.8 s at maxm = 30 / 60 / 240
  vev, T=0.1 (beta=10):   8.5 s / 12.4 s at maxm = 30 / 60, and the answer
                          moves by 6e-11 between them

i.e. the real-time correlator saturates whatever cap it is given, while the
static one saturates its own entanglement below chi=60 -- METTS states are
minimally entangled by construction, which is the method's selling point
and also why `vev` cannot reach this port's chi ~ 120-160 crossover on any
hardware. That is the deliverable for `vev`: a bound, not a speedup.

The bond dimension a METTS run actually reaches is also unusual in a way
that matters for the dispatch knobs: every sample restarts from a fresh
bond-dimension-1 product state and grows, so one run walks the *whole*
ladder of bond dimensions rather than sitting at one -- a shape zoo that
XLA compiles a kernel per rung of (see vevtk/mettsvev.py's own note on why
the old per-call numba/JAX kernels made METTS slower, which is the same
observation from the other side). backend.set_pad_bonds is the direct
answer to that, so --pad-bonds is expected to matter *more* here than for
a calculation that stays at one bond dimension, and both passes are run.

Same conventions as the rest of this directory: fixed seed, cold and warm
reported separately (XLA compiles per shape, so a one-shot script and a
long-running process see very different costs), agreement against the
host result reported as a column rather than assumed, and the CPU
baseline submitted as its own job on a CPU partition (CPU-only work does
not belong on a GPU partition), merged afterwards with --from-json.

METTS is stochastic, but not *here*: a fixed seed fixes the whole Markov
trajectory, so host and device must return the same number to ~1e-14
unless a 1e-15 perturbation flips an rng.choice draw -- in which case the
deviation column blows up rather than drifting, which is the failure this
column is for.
"""

import argparse
import json
import os
import sys
import time


def _exchange(sc, i, j):
    return sc.Sx[i] * sc.Sx[j] + sc.Sy[i] * sc.Sy[j] + sc.Sz[i] * sc.Sz[j]


def build_chain(n, maxm, seed, field=0.0):
    """A uniform S=1/2 Heisenberg chain, the same model the rest of this
    directory uses. A small uniform field is available so the static
    observable has a non-trivial thermal value to agree on."""
    import numpy as np
    from dmrgpy import spinchain

    np.random.seed(seed)
    sc = spinchain.Spin_Chain(["S=1/2" for _ in range(n)],
                              itensor_version="python")
    h = 0
    for i in range(n - 1):
        h = h + _exchange(sc, i, i + 1)
    if field:
        for i in range(n):
            h = h + field * sc.Sz[i]
    sc.set_hamiltonian(h)
    sc.maxm = maxm
    return sc


def time_metts_vev(args, maxm):
    """(times, value) for `reps` runs of one metts_vev sampling run."""
    times, value = [], None
    for _ in range(args.reps):
        sc = build_chain(args.n, maxm, args.seed, field=args.field)
        t0 = time.perf_counter()
        mean, _err = sc.metts_vev(sc.Sz[args.n // 2], T=args.T,
                                  nsamples=args.nsamples,
                                  nwarmup=args.nwarmup,
                                  dbeta_half_step=args.dbeta,
                                  seed=args.seed)
        times.append(time.perf_counter() - t0)
        value = complex(mean)
    return times, value


def time_metts_dc(args, maxm):
    """(times, C(t)) for `reps` runs of one dynamical-METTS correlator."""
    import numpy as np
    times, value = [], None
    for _ in range(args.reps):
        sc = build_chain(args.n, maxm, args.seed, field=args.field)
        t0 = time.perf_counter()
        ts, cs, _errs = sc.metts_dynamical_correlator(
                name=(sc.Sz[args.n // 2], sc.Sz[args.n // 2]), T=args.T,
                nt=args.nt, dt=args.dt, nsamples=args.nsamples,
                nwarmup=args.nwarmup, dbeta_half_step=args.dbeta,
                seed=args.seed)
        times.append(time.perf_counter() - t0)
        value = np.asarray(cs)
    return times, value


CALCS = {"vev": time_metts_vev, "dc": time_metts_dc}


def _deviation(value, reference):
    import numpy as np
    if reference is None:
        return 0.0
    return float(np.max(np.abs(np.asarray(value) - np.asarray(reference))))


def run_backend(name, args, results):
    from dmrgpy.pyitensor import backend as bk

    try:
        bk.set_backend(name)
    except Exception as exc:
        print("[%s] unavailable: %s" % (name, exc))
        return
    print("\n=== backend %s: %s ===" % (name, bk.device_info()))
    results["meta"]["device_" + name] = bk.device_info()
    results["metts"][name] = {}

    for calc in args.calcs.split(","):
        if calc not in CALCS:
            raise SystemExit("unknown --calcs entry %r (have: %s)"
                             % (calc, ", ".join(CALCS)))
        reference = None
        for maxm in [int(m) for m in args.maxm.split(",")]:
            if args.pad_bonds:
                bk.set_pad_bonds(maxm)
            try:
                times, value = CALCS[calc](args, maxm)
            except Exception as exc:
                print("  %-4s maxm=%-4d FAILED %s: %s"
                      % (calc, maxm, type(exc).__name__, exc))
                results["metts"][name].setdefault(calc, {})[str(maxm)] = {
                        "error": "%s: %s" % (type(exc).__name__, exc)}
                bk.set_pad_bonds(None)
                continue
            warm = min(times[1:]) if len(times) > 1 else times[0]
            # The reference is this backend's own smallest-maxm result: a
            # cross-backend comparison is the report's job (it needs both
            # JSONs), while what a single job can check on its own is that
            # nothing in the sweep returned a wildly different trajectory.
            dev = _deviation(value, reference)
            if reference is None:
                reference = value
            entry = {"times": times, "cold": times[0], "warm": warm,
                     "max_dev_vs_first": dev,
                     "compilations": bk.compilations(),
                     "value": (str(value) if calc == "vev"
                               else [str(v) for v in value])}
            results["metts"][name].setdefault(calc, {})[str(maxm)] = entry
            print("  %-4s maxm=%-4d cold %8.2fs  warm %8.2fs  dev %.2e  "
                  "kernels %s" % (calc, maxm, times[0], warm, dev,
                                  entry["compilations"]))
            sys.stdout.flush()
            bk.set_pad_bonds(None)


def report(results):
    print("\n" + "=" * 78)
    meta = results.get("meta", {})
    print("METTS: n=%s, T=%s, nsamples=%s, nwarmup=%s, dbeta=%s, nt=%s, dt=%s"
          % (meta.get("n"), meta.get("T"), meta.get("nsamples"),
             meta.get("nwarmup"), meta.get("dbeta"), meta.get("nt"),
             meta.get("dt")))
    print("=" * 78)
    block = results.get("metts", {})
    calcs = sorted({c for b in block.values() for c in b})
    for calc in calcs:
        maxms = sorted({m for b in block.values() for m in b.get(calc, {})},
                       key=int)
        print("\n%s" % {"vev": "metts_vev (imaginary time only)",
                        "dc": "metts_dynamical_correlator (+ real time)"}[calc])
        print("  %-12s %-8s %12s %12s %10s" % ("backend", "maxm", "cold(s)",
                                               "warm(s)", "kernels"))
        for name in sorted(block):
            for maxm in maxms:
                e = block[name].get(calc, {}).get(maxm)
                if not e:
                    continue
                if "error" in e:
                    print("  %-12s %-8s %s" % (name, maxm, e["error"][:48]))
                    continue
                print("  %-12s %-8s %12.2f %12.2f %10s"
                      % (name, maxm, e["cold"], e["warm"],
                         e["compilations"]))
        host_keys = sorted(k for k in block if k.startswith("numpy"))
        dev_keys = sorted(k for k in block if k.startswith("jax"))
        for maxm in maxms:
            for hk in host_keys:
                h = block[hk].get(calc, {}).get(maxm)
                if not h or "error" in h:
                    continue
                for dk in dev_keys:
                    d = block[dk].get(calc, {}).get(maxm)
                    if not d or "error" in d:
                        continue
                    print("  -> maxm=%-5s %s vs %s: %.2fx cold, %.2fx warm"
                          % (maxm, dk, hk, h["cold"] / d["cold"],
                             h["warm"] / d["warm"]))
        # Cross-backend agreement, the column that says the two runs
        # followed the same Markov trajectory rather than merely landing
        # near each other.
        import numpy as np
        for maxm in maxms:
            vals = {}
            for name in block:
                e = block[name].get(calc, {}).get(maxm)
                if e and "error" not in e:
                    v = e["value"]
                    vals[name] = (complex(v) if isinstance(v, str)
                                  else np.array([complex(x) for x in v]))
            if len(vals) > 1:
                ref_name = sorted(vals)[0]
                ref = vals[ref_name]
                for name, v in sorted(vals.items()):
                    if name == ref_name:
                        continue
                    print("  -> maxm=%-5s %s vs %s: max|dC| = %.2e"
                          % (maxm, name, ref_name,
                             float(np.max(np.abs(v - ref)))))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--src", default=None,
                    help="checkout's src/ to put at the front of sys.path")
    ap.add_argument("--backends", default="numpy",
                    help="comma-separated: numpy, jax")
    ap.add_argument("--calcs", default="vev,dc",
                    help="comma-separated: vev (metts_vev), dc "
                         "(metts_dynamical_correlator)")
    ap.add_argument("--n", type=int, default=20, help="chain length")
    ap.add_argument("--maxm", default="30,60,120,240",
                    help="MPS bond dimensions to sweep")
    ap.add_argument("--T", type=float, default=0.1,
                    help="temperature (beta=1/T). Default beta=10: a real "
                         "low-temperature METTS run, i.e. the most entangled "
                         "the sampled states realistically get -- see the "
                         "module docstring on why T=0.5 makes the sweep flat")
    ap.add_argument("--field", type=float, default=0.3,
                    help="uniform Sz field, so <Sz> has a non-trivial value")
    ap.add_argument("--nsamples", type=int, default=2,
                    help="retained METTS samples -- kept small on purpose: "
                         "cost is linear in it and this measures per-sample "
                         "cost, not statistical quality")
    ap.add_argument("--nwarmup", type=int, default=1)
    ap.add_argument("--dbeta", type=float, default=0.1,
                    help="imaginary-time step for the e^{-beta H/2} evolution")
    ap.add_argument("--nt", type=int, default=100,
                    help="real-time steps (dc). Default 100, i.e. t <= 10 at "
                         "dt=0.1: shorter than that and the evolved states "
                         "never reach the bond-dimension cap")
    ap.add_argument("--dt", type=float, default=0.1)
    ap.add_argument("--reps", type=int, default=2,
                    help=">=2 so the cold (compiling) run and the warm runs "
                         "are separated")
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument("--pad-bonds", action="store_true",
                    help="freeze every MPS bond at the sweep's bond dimension "
                         "(backend.set_pad_bonds), which is also what turns "
                         "backend.set_jit(\"auto\") on -- expected to matter "
                         "more for METTS than for anything else here, since "
                         "every sample regrows chi from 1")
    ap.add_argument("--json", default=None)
    ap.add_argument("--from-json", default=None, nargs="*",
                    help="merge previous runs' JSON and re-print the report "
                         "without recomputing (how the CPU and GPU jobs are "
                         "combined)")
    args = ap.parse_args()

    if args.src:
        sys.path.insert(0, os.path.abspath(args.src))

    if args.from_json:
        # Same rule as tdz_bench.py: a padded pass and an eager pass are
        # both "jax" and share maxm keys, so merging on the bare backend
        # name would have one silently replace the other.
        merged = {"meta": {}, "metts": {}}
        for path in args.from_json:
            with open(path) as fh:
                part = json.load(fh)
            suffix = "+pad" if part.get("meta", {}).get("pad_bonds") else ""
            for name, calcs in part.get("metts", {}).items():
                dest = merged["metts"].setdefault(name + suffix, {})
                for calc, rows in calcs.items():
                    dest.setdefault(calc, {}).update(rows)
            for mk, mv in part.get("meta", {}).items():
                if mk in ("pad_bonds", "maxm", "device_jax", "device_numpy"):
                    merged["meta"].setdefault("per_run", {}).setdefault(
                            os.path.basename(path), {})[mk] = mv
                else:
                    merged["meta"].setdefault(mk, mv)
        report(merged)
        return

    results = {"meta": {k: getattr(args, k) for k in
                        ("n", "maxm", "T", "field", "nsamples", "nwarmup",
                         "dbeta", "nt", "dt", "reps", "seed", "pad_bonds",
                         "calcs")},
               "metts": {}}
    for name in args.backends.split(","):
        run_backend(name.strip(), args, results)
    report(results)
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(results, fh, indent=1)
        print("\nwrote %s" % args.json)


if __name__ == "__main__":
    main()
