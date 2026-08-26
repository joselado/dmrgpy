"""Phase 0 of the pyitensor GPU port (docs/pyitensor_gpu_port_plan.md):
run *unmodified* pyitensor on a GPU node and see what actually happens.

Two questions, neither of which the primitive-level microbenchmark
(gpu_microbench.py, run alongside this) can answer:

1. pyitensor/kernels.py's USE_JAX defaults to _detect_default_use_jax(),
   which auto-enables the fused-einsum JAX matvec whenever jax reports any
   non-CPU device -- a heuristic its own docstring calls UNTESTED on real
   GPU hardware. On the cluster it fires, because the cluster's default python3
   ships jax with its CUDA plugin. So every user running pyitensor on a GPU
   node today silently takes that path. This script forces USE_JAX both
   ways over the same chains and reports the ratio, which either vindicates
   the heuristic or condemns it.
2. Whether the JAX path is transfer-bound as predicted: it converts
   numpy<->jax per matvec call, so on a GPU each call pays an H2D+D2H round
   trip. Compare the per-call transfer cost the microbenchmark measures
   against its matvec cost at the same chi.

Ground-state energies from both paths are compared against each other
(they must agree to DMRG's own tolerance -- same algorithm, different
contraction backend) so a timing win can't be a silently wrong answer.
"""

import argparse
import json
import os
import sys
import time


def build_chain(n, maxm, nsweeps):
    """Uniform S=1/2 Heisenberg chain on the pure-Python backend -- the
    same model benchmarks/run_benchmarks.py sweeps, so numbers here are
    comparable to that harness."""
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
    sc.nsweeps = nsweeps
    return sc


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--src", default=None,
                    help="path to dmrgpy's src/ (prepended to sys.path)")
    ap.add_argument("--sizes", default="16:60,20:100,24:200",
                    help="comma-separated n:maxm pairs")
    ap.add_argument("--nsweeps", type=int, default=4)
    ap.add_argument("--paths", default="plain,jax",
                    help="which matvec paths to time: plain, jax, or both "
                         "(the CPU baseline job wants plain only -- the JAX "
                         "path's CPU regression is already known)")
    ap.add_argument("--json", default=None)
    args = ap.parse_args()

    if args.src:
        sys.path.insert(0, os.path.abspath(args.src))

    out = {"meta": {}, "runs": []}

    # Report the environment before touching dmrgpy, so a failure below
    # still leaves the diagnostic behind.
    try:
        import jax
        jax.config.update("jax_enable_x64", True)
        out["meta"]["jax_version"] = jax.__version__
        out["meta"]["jax_devices"] = [str(d) for d in jax.devices()]
        out["meta"]["jax_default_platform"] = jax.devices()[0].platform
    except Exception as exc:
        out["meta"]["jax_error"] = "%s: %s" % (type(exc).__name__, exc)
    try:
        import torch
        out["meta"]["torch_version"] = torch.__version__
        out["meta"]["torch_cuda"] = bool(torch.cuda.is_available())
        if torch.cuda.is_available():
            out["meta"]["torch_device"] = torch.cuda.get_device_name(0)
    except Exception as exc:
        out["meta"]["torch_error"] = "%s: %s" % (type(exc).__name__, exc)

    from dmrgpy.pyitensor import kernels
    out["meta"]["dmrgpy_path"] = os.path.dirname(
        sys.modules["dmrgpy"].__file__)
    out["meta"]["kernels_use_jax_default"] = bool(kernels.USE_JAX)
    out["meta"]["kernels_available"] = kernels.available() \
        if hasattr(kernels, "available") else None
    print("jax:   %s" % out["meta"].get("jax_devices"))
    print("torch: %s" % out["meta"].get("torch_device",
                                        out["meta"].get("torch_cuda")))
    print("kernels.USE_JAX default on this node: %s"
          % out["meta"]["kernels_use_jax_default"])
    print("(that default is the hazard this script exists to test)\n")

    for spec in args.sizes.split(","):
        n, maxm = (int(x) for x in spec.split(":"))
        row = {"n": n, "maxm": maxm, "nsweeps": args.nsweeps}
        paths = [p.strip() for p in args.paths.split(",") if p.strip()]
        for use_jax in [p == "jax" for p in paths]:
            kernels.USE_JAX = use_jax
            tag = "jax" if use_jax else "plain"
            try:
                sc = build_chain(n, maxm, args.nsweeps)
                t0 = time.perf_counter()
                e = sc.gs_energy(mode="DMRG")
                dt = time.perf_counter() - t0
                row[tag] = {"seconds": dt, "energy": float(e.real)
                            if hasattr(e, "real") else float(e)}
                print("  n=%-3d maxm=%-4d %-5s  %8.2f s   E=%.10f"
                      % (n, maxm, tag, dt, row[tag]["energy"]))
            except Exception as exc:
                row[tag] = {"error": "%s: %s" % (type(exc).__name__, exc)}
                print("  n=%-3d maxm=%-4d %-5s  FAILED %s"
                      % (n, maxm, tag, row[tag]["error"]))
            sys.stdout.flush()

        if "seconds" in row.get("plain", {}) and "seconds" in row.get("jax", {}):
            row["jax_speedup"] = row["plain"]["seconds"] / row["jax"]["seconds"]
            row["energy_diff"] = abs(row["plain"]["energy"]
                                     - row["jax"]["energy"])
            print("    -> jax speedup %.2fx, |dE| = %.2e"
                  % (row["jax_speedup"], row["energy_diff"]))
        out["runs"].append(row)
        sys.stdout.flush()

    if args.json:
        with open(args.json, "w") as f:
            json.dump(out, f, indent=1)
        print("\nwrote %s" % args.json)


if __name__ == "__main__":
    main()
