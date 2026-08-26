"""Turn the Phase 0 JSON files (gpu_microbench.py / pyitensor_gpu_probe.py,
run by phase0_gpu.the scheduler and phase0_cpu.the scheduler) into the two things the
decision in docs/pyitensor_gpu_port_plan.md actually needs: a plot of
where the GPU crosses over, and a printed verdict against the plan's own
gate (no crossover below chi ~ 512 on an A100 -> report, don't port).

Plots rather than only printing, per CLAUDE.md's rule for anything in this
repo that produces a sequence of values worth looking at. Three panels:

1. absolute time vs chi, per op, CPU and GPU overlaid (log-log) -- the raw
   measurement;
2. speedup vs chi, with the break-even line at 1.0 -- the crossover;
3. matvec time vs transfer time on the GPU -- whether a per-call
   conversion design (what pyitensor/kernels.py's JAX path does today)
   could ever win, independently of how fast the device is.

Usage:
    python3 phase0_report.py --cpu micro-cpu1-*.json --gpu micro-gpu-*.json
                             [--probe probe-*.json] [--out phase0.png]
"""

import argparse
import glob
import json
import os


def load(patterns):
    """Read every JSON matching any pattern, newest last (so a rerun wins)."""
    paths = []
    for pat in patterns:
        paths.extend(sorted(glob.glob(pat), key=os.path.getmtime))
    out = []
    for p in paths:
        with open(p) as f:
            out.append((p, json.load(f)))
    return out


def series(doc, backend, op):
    """{chi: median seconds} for one backend/op, skipping failed points."""
    data = doc.get("data", {}).get(backend, {})
    ops = data.get("ops", {})
    got = {}
    for chi, rec in ops.get(op, {}).items():
        if isinstance(rec, dict) and "median_s" in rec:
            got[int(chi)] = rec["median_s"]
    return got


def backends_in(doc):
    return [b for b, v in doc.get("data", {}).items()
            if isinstance(v, dict) and v.get("ops")]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cpu", nargs="*", default=[],
                    help="microbench JSON(s) from the CPU job")
    ap.add_argument("--gpu", nargs="*", default=[],
                    help="microbench JSON(s) from the GPU job")
    ap.add_argument("--probe", nargs="*", default=[],
                    help="pyitensor_gpu_probe JSON(s)")
    ap.add_argument("--gate-chi", type=int, default=512,
                    help="the plan's gate: a crossover at or below this chi")
    ap.add_argument("--out", default="phase0.png")
    args = ap.parse_args()

    cpu_docs, gpu_docs = load(args.cpu), load(args.gpu)
    if not cpu_docs or not gpu_docs:
        raise SystemExit("need at least one CPU and one GPU microbench JSON")

    cpu = cpu_docs[-1][1]
    gpu = gpu_docs[-1][1]
    gpu_backends = backends_in(gpu)
    ops = ["matvec", "eigh", "svd", "qr"]

    print("CPU: %s" % cpu["meta"].get("platform"))
    print("GPU: %s" % (gpu["meta"].get("nvidia_smi") or "?"))
    print("GPU backends measured: %s\n" % ", ".join(gpu_backends))

    verdict = {}
    for gb in gpu_backends:
        print("=== %s vs numpy (CPU) ===" % gb)
        print("%-8s %8s %12s %12s %9s" % ("op", "chi", "cpu ms", "gpu ms", "speedup"))
        for op in ops:
            c, g = series(cpu, "numpy", op), series(gpu, gb, op)
            crossover = None
            for chi in sorted(set(c) & set(g)):
                sp = c[chi] / g[chi]
                print("%-8s %8d %12.3f %12.3f %8.2fx"
                      % (op, chi, 1e3 * c[chi], 1e3 * g[chi], sp))
                if sp > 1.0 and crossover is None:
                    crossover = chi
            verdict[(gb, op)] = crossover
            print("  -> %s crossover: %s\n"
                  % (op, "chi=%d" % crossover if crossover else "none in range"))

        # The transfer question: is the matvec even worth shipping?
        mv, tr = series(gpu, gb, "matvec"), series(gpu, gb, "transfer")
        both = sorted(set(mv) & set(tr))
        if both:
            print("  per-call transfer vs matvec on the device (%s):" % gb)
            for chi in both:
                print("    chi=%-5d matvec %8.3f ms   round trip %8.3f ms   %s"
                      % (chi, 1e3 * mv[chi], 1e3 * tr[chi],
                         "transfer-bound" if tr[chi] > mv[chi] else "compute-bound"))
            print()

    for path, doc in load(args.probe):
        print("=== probe: %s ===" % os.path.basename(path))
        print("  jax devices: %s" % doc["meta"].get("jax_devices"))
        print("  kernels.USE_JAX default here: %s"
              % doc["meta"].get("kernels_use_jax_default"))
        for row in doc.get("runs", []):
            bits = ["n=%d maxm=%d" % (row["n"], row["maxm"])]
            for tag in ("plain", "jax"):
                r = row.get(tag)
                if not r:
                    continue
                bits.append("%s=%s" % (tag, "%.2fs" % r["seconds"]
                                       if "seconds" in r else r.get("error")))
            if "jax_speedup" in row:
                bits.append("jax %.2fx, |dE|=%.1e"
                            % (row["jax_speedup"], row["energy_diff"]))
            print("  " + "  ".join(bits))
        print()

    # ---- the gate ----
    best = [chi for (gb, op), chi in verdict.items()
            if op == "matvec" and chi is not None]
    print("=" * 68)
    if best and min(best) <= args.gate_chi:
        print("GATE PASSED: matvec crosses over at chi=%d (<= %d). Proceed to "
              "Phase 1." % (min(best), args.gate_chi))
    elif best:
        print("GATE MARGINAL: matvec crosses over only at chi=%d (> %d). "
              "Report before porting." % (min(best), args.gate_chi))
    else:
        print("GATE FAILED: no matvec crossover in the measured range. "
              "Report, do not port.")
    print("=" * 68)

    # ---- plot ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("\nmatplotlib unavailable; skipped the figure")
        return

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))

    ax = axes[0]
    for op, style in zip(ops, ["o-", "s-", "^-", "v-"]):
        c = series(cpu, "numpy", op)
        if c:
            ax.plot(sorted(c), [1e3 * c[k] for k in sorted(c)], style,
                    color="0.35", label="CPU %s" % op)
        for gb in gpu_backends:
            g = series(gpu, gb, op)
            if g:
                ax.plot(sorted(g), [1e3 * g[k] for k in sorted(g)], style,
                        label="%s %s" % (gb, op))
    ax.set_xscale("log", base=2); ax.set_yscale("log")
    ax.set_xlabel(r"bond dimension $\chi$"); ax.set_ylabel("time per call [ms]")
    ax.set_title("primitive cost, complex128")
    ax.legend(fontsize=6, ncol=2); ax.grid(alpha=.3)

    ax = axes[1]
    for gb in gpu_backends:
        for op, style in zip(ops, ["o-", "s-", "^-", "v-"]):
            c, g = series(cpu, "numpy", op), series(gpu, gb, op)
            chis = sorted(set(c) & set(g))
            if chis:
                ax.plot(chis, [c[k] / g[k] for k in chis], style,
                        label="%s %s" % (gb, op))
    ax.axhline(1.0, color="k", lw=1, ls="--")
    ax.set_xscale("log", base=2); ax.set_yscale("log")
    ax.set_xlabel(r"bond dimension $\chi$"); ax.set_ylabel("GPU speedup over CPU")
    ax.set_title("crossover (above the dashed line the GPU wins)")
    ax.legend(fontsize=6, ncol=2); ax.grid(alpha=.3)

    ax = axes[2]
    for gb in gpu_backends:
        mv, tr = series(gpu, gb, "matvec"), series(gpu, gb, "transfer")
        chis = sorted(set(mv) & set(tr))
        if chis:
            ax.plot(chis, [1e3 * mv[k] for k in chis], "o-",
                    label="%s matvec" % gb)
            ax.plot(chis, [1e3 * tr[k] for k in chis], "s--",
                    label="%s H2D+D2H" % gb)
    ax.set_xscale("log", base=2); ax.set_yscale("log")
    ax.set_xlabel(r"bond dimension $\chi$"); ax.set_ylabel("time [ms]")
    ax.set_title("why arrays must stay resident")
    ax.legend(fontsize=7); ax.grid(alpha=.3)

    fig.tight_layout()
    fig.savefig(args.out, dpi=140)
    print("\nwrote %s" % args.out)


if __name__ == "__main__":
    main()
