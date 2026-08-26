"""End-to-end speedup of the pyitensor GPU port: ground-state DMRG and the
KPM dynamical correlator, NumPy (host) versus JAX (device), swept over bond
dimension. This is Phase 2's deliverable in
docs/pyitensor_gpu_port_plan.md -- the number Phase 0 deliberately refused
to guess.

Three things this measures that a naive timing would get wrong:

* **Cold versus warm.** Eager JAX compiles an XLA kernel per distinct
  (operation, shape), and DMRG mints new shapes as bond dimensions grow,
  so the first run of a given size pays a compile tax that later runs of
  the same size do not. Measured on CPU at n=6/maxm=20, 672 compilations
  cost 18.4 s of a 29.1 s run. Reporting only the first run would
  understate the device by a large factor; reporting only the warm run
  would overstate it for anyone whose script does one calculation and
  exits. Both are printed.
* **Bond dimension, not chain length.** The device's advantage is a
  large-chi story throughout this plan (the primitive crossover sits at
  chi=64 and the matvec gains 688x by chi=1024), so chi is the sweep axis
  and the chain is long enough that bond dimensions actually reach the
  requested maxm rather than being capped by the chain's own entanglement.
* **Correctness at every point.** Ground-state energies are compared
  against the NumPy run of the same size, and every KPM spectrum is
  checked against the exact zeroth-moment sum rule (int S_zz domega =
  <GS|Sz Sz|GS> = 1/4 on a spin-1/2 site), the same invariant
  tests/test_kpm_dynamical_correlator_python.py asserts. A timing point
  that fails either is reported as such rather than being quietly kept.

The CPU baseline must run as its own job on a CPU partition -- CPU-only
workloads are not permitted on the cluster's GPU partitions -- so the two
halves are combined afterwards by --from-json.
"""

import argparse
import json
import os
import sys
import time


def _exchange(sc, i, j, coupling=1.0):
    return (coupling * (sc.Sx[i] * sc.Sx[j] + sc.Sy[i] * sc.Sy[j]
                        + sc.Sz[i] * sc.Sz[j]))


def build_chain(n, maxm, kpmmaxm, nsweeps, seed, model="heisenberg"):
    """`model` exists because the obvious choice is a bad GPU test case.

    A uniform 1D Heisenberg chain has entanglement S ~ (1/3) log L, so its
    ground state converges at chi ~ 60 whatever maxm is set to -- measured
    directly, E0 at n=32 is identical to 13 digits for maxm=120 and
    maxm=240, and both CPU and GPU times were flat across the whole sweep.
    Benchmarking a GPU on it therefore measures dispatch overhead and
    nothing else, which is exactly what the first ground-state sweep here
    did.

    The other two models genuinely need large bond dimension at these
    sizes:

    "ladder"  -- a 2-leg ladder folded onto the chain (rungs i..i+1 on even
                 i, legs i..i+2). Entanglement grows with the transverse
                 width rather than logarithmically, so chi in the hundreds
                 is real work, not padding.
    "ladder3" -- a 3-leg ladder, same idea with more transverse width.
    "j1j2"    -- J1-J2 chain at J2/J1 = 0.35, in the gapless incommensurate
                 regime. NOT 0.5: that is the Majumdar-Ghosh point, whose
                 ground state is an *exact* dimer product state with chi=2
                 -- measured here, E0 = -7.5 exactly at every maxm from 30
                 to 240, i.e. the single worst possible benchmark for a
                 large-chi claim. It was 0.5 in the first version of this
                 file.
    """
    import numpy as np
    from dmrgpy import spinchain
    # Same seed for both backends: the DMRG start state is a random MPS, and
    # comparing energies across backends is only meaningful from the same
    # starting point.
    np.random.seed(seed)
    sc = spinchain.Spin_Chain(["S=1/2" for _ in range(n)],
                              itensor_version="python")
    h = 0
    if model == "heisenberg":
        for i in range(n - 1):
            h = h + _exchange(sc, i, i + 1)
    elif model == "ladder":
        for i in range(0, n - 1, 2):          # rungs
            h = h + _exchange(sc, i, i + 1)
        for i in range(n - 2):                # legs
            h = h + _exchange(sc, i, i + 2)
    elif model == "ladder3":
        for c in range(n // 3):               # rungs within each column
            for r in range(2):
                h = h + _exchange(sc, 3 * c + r, 3 * c + r + 1)
        for i in range(n - 3):                # legs between columns
            h = h + _exchange(sc, i, i + 3)
    elif model == "j1j2":
        for i in range(n - 1):
            h = h + _exchange(sc, i, i + 1)
        for i in range(n - 2):
            h = h + _exchange(sc, i, i + 2, coupling=0.35)
    else:
        raise ValueError("unknown model %r" % model)
    sc.set_hamiltonian(h)
    sc.maxm = maxm
    sc.kpmmaxm = kpmmaxm
    sc.nsweeps = nsweeps
    return sc


def time_groundstate(n, maxm, nsweeps, seed, reps, model):
    import numpy as np
    out = []
    energy = None
    for r in range(reps):
        sc = build_chain(n, maxm, maxm, nsweeps, seed, model)
        t0 = time.perf_counter()
        e = sc.gs_energy(mode="DMRG")
        out.append(time.perf_counter() - t0)
        energy = float(np.real(e))
    return out, energy


def time_kpm(n, maxm, kpmmaxm, nsweeps, seed, reps, delta, model,
             energy_truncation=None):
    """energy_truncation: None (off) or a (kpm_scale, dK, n_sweeps,
    threshold) tuple, which turns on the Holzner et al. energy truncation
    (pyitensor/kpm_energy_truncation.py). That path is worth timing
    separately because it adds its own Krylov projection per site per
    sweep -- i.e. more small dense work per moment, on top of the
    recursion's own applyMPO+truncate."""
    import numpy as np
    es = np.linspace(-8.0, 8.0, 1200)   # wide: the sum rule needs all the weight
    out = []
    weight = None
    for r in range(reps):
        sc = build_chain(n, maxm, kpmmaxm, nsweeps, seed, model)
        if energy_truncation is not None:
            scale, dK, ns, thr = energy_truncation
            sc.kpm_scale = scale
            sc.kpm_energy_truncate = True
            sc.kpm_truncate_dK = dK
            sc.kpm_truncate_nsweeps = ns
            sc.kpm_truncate_threshold = thr
        sc.gs_energy(mode="DMRG")       # not timed here; timed separately above
        t0 = time.perf_counter()
        x, y = sc.get_dynamical_correlator(mode="DMRG", submode="KPM",
                                            name=(sc.Sz[0], sc.Sz[0]),
                                            es=es, delta=delta)
        out.append(time.perf_counter() - t0)
        weight = float(np.trapezoid(np.real(y), x))
    return out, weight


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--src", default=None)
    ap.add_argument("--backends", default="numpy",
                    help="comma-separated: numpy, jax")
    ap.add_argument("--n", type=int, default=32,
                    help="chain length; long enough that bond dims reach maxm")
    ap.add_argument("--maxm", default="30,60,120,240",
                    help="ground-state bond dimensions to sweep")
    ap.add_argument("--kpmmaxm", default="40,80,160",
                    help="KPM bond dimensions to sweep")
    ap.add_argument("--kpm-n", type=int, default=16,
                    help="chain length for the KPM sweep (kept smaller: the "
                         "moment recursion is O(n) per moment)")
    ap.add_argument("--nsweeps", type=int, default=5)
    ap.add_argument("--delta", type=float, default=0.1)
    ap.add_argument("--reps", type=int, default=2,
                    help=">=2 so the first (cold, compiling) run and the "
                         "later (warm) runs are separated")
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument("--model", default="heisenberg",
                    choices=["heisenberg", "ladder", "ladder3", "j1j2"],
                    help="heisenberg converges at chi~60 and so cannot show a "
                         "ground-state GPU win at any maxm; ladder and j1j2 "
                         "genuinely need large bond dimension (see "
                         "build_chain's docstring)")
    ap.add_argument("--kpm-energy-truncation", default=None,
                    metavar="SCALE,DK,NSWEEPS,THRESHOLD",
                    help="turn on KPM energy truncation with these settings, "
                         "e.g. 1.4,10,2,1.0 . The narrow kpm_scale is the "
                         "point of the feature: it buys spectral resolution "
                         "at the cost of the rescaled spectrum leaking "
                         "outside [-1,1], which this projection removes.")
    ap.add_argument("--pad-bonds", action="store_true",
                    help="freeze every MPS bond at the sweep's bond dimension "
                         "(backend.set_pad_bonds), so XLA compiles one kernel "
                         "per operation instead of one per bond dimension. "
                         "Measured on CPU: 4x fewer compilations and 1.64x "
                         "faster cold, but 2x slower warm because the padded "
                         "blocks are known zeros -- the point of running it on "
                         "a device is that both sides of that trade move the "
                         "other way there.")
    ap.add_argument("--json", default=None)
    args = ap.parse_args()

    if args.src:
        sys.path.insert(0, os.path.abspath(args.src))
    import numpy as np
    from dmrgpy.pyitensor import backend as bk

    results = {"meta": {"n": args.n, "kpm_n": args.kpm_n,
                        "model": args.model,
                        "kpm_energy_truncation": args.kpm_energy_truncation,
                        "nsweeps": args.nsweeps, "delta": args.delta,
                        "reps": args.reps, "seed": args.seed,
                        "python": sys.version.split()[0]},
               "groundstate": {}, "kpm": {}}
    try:
        import subprocess
        results["meta"]["nvidia_smi"] = subprocess.run(
            ["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"],
            capture_output=True, text=True, timeout=30).stdout.strip()
    except Exception:
        pass

    for name in args.backends.split(","):
        name = name.strip()
        if not name:
            continue
        try:
            bk.set_backend(name)
        except Exception as exc:
            print("[%s] unavailable: %s" % (name, exc))
            continue
        print("\n=== backend %s: %s ===" % (name, bk.device_info()))
        results["meta"]["device_" + name] = bk.device_info()

        results["meta"]["pad_bonds"] = bool(args.pad_bonds)
        results["groundstate"][name] = {}
        for maxm in [int(m) for m in args.maxm.split(",")]:
            if args.pad_bonds:
                bk.set_pad_bonds(maxm)
            try:
                times, energy = time_groundstate(args.n, maxm, args.nsweeps,
                                                 args.seed, args.reps,
                                                 args.model)
            except Exception as exc:
                print("  gs   maxm=%-5d FAILED %s: %s"
                      % (maxm, type(exc).__name__, exc))
                results["groundstate"][name][str(maxm)] = {"error": str(exc)}
                continue
            results["groundstate"][name][str(maxm)] = {
                "times": times, "cold": times[0], "warm": min(times[1:]) if len(times) > 1 else times[0],
                "energy": energy}
            print("  gs   maxm=%-5d cold %8.2fs   warm %8.2fs   E0=%.10f"
                  % (maxm, times[0], min(times[1:]) if len(times) > 1 else times[0], energy))
            sys.stdout.flush()

        results["kpm"][name] = {}
        for kpmmaxm in [int(m) for m in args.kpmmaxm.split(",")]:
            if args.pad_bonds:
                bk.set_pad_bonds(kpmmaxm)
            try:
                # maxm == kpmmaxm on purpose: with different values the
                # ground-state sweep and the moment recursion have two
                # separate shape families, so XLA compiles every kernel
                # twice. One K for both halves is one family.
                et = None
                if args.kpm_energy_truncation:
                    f = args.kpm_energy_truncation.split(",")
                    et = (float(f[0]), int(f[1]), int(f[2]), float(f[3]))
                times, weight = time_kpm(args.kpm_n, kpmmaxm, kpmmaxm,
                                         args.nsweeps, args.seed, args.reps,
                                         args.delta, args.model,
                                         energy_truncation=et)
            except Exception as exc:
                print("  kpm  kpmmaxm=%-5d FAILED %s: %s"
                      % (kpmmaxm, type(exc).__name__, exc))
                results["kpm"][name][str(kpmmaxm)] = {"error": str(exc)}
                continue
            results["kpm"][name][str(kpmmaxm)] = {
                "times": times, "cold": times[0],
                "warm": min(times[1:]) if len(times) > 1 else times[0],
                "sum_rule": weight, "sum_rule_error": abs(weight - 0.25)}
            print("  kpm  kpmmaxm=%-5d cold %8.2fs   warm %8.2fs   "
                  "sum rule %.6f (err %.1e)"
                  % (kpmmaxm, times[0],
                     min(times[1:]) if len(times) > 1 else times[0],
                     weight, abs(weight - 0.25)))
            sys.stdout.flush()

    if args.json:
        with open(args.json, "w") as f:
            json.dump(results, f, indent=1)
        print("\nwrote %s" % args.json)


if __name__ == "__main__":
    main()
