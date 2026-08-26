"""What `backend.set_jit` buys: the dispatch floor, measured.

This is item 1 of docs/pyitensor_gpu_port_plan.md Sec. 9 -- the one item
on that list that is not a port. Nothing here changes what the engine
computes; it changes how many kernels the engine *dispatches* to get
there, which is the quantity that decides whether a device can win at all
below chi ~ 120-160.

The two knobs only work together, and this script is built to show that
rather than assert it:

* eager, unpadded -- the port's original behavior: three to twelve
  dispatches per matvec, five per contraction, one per Krylov basis
  vector, each paying a floor of ~0.35 ms on an H200 whatever its size;
* padded, eager -- one shape per operation, so XLA stops compiling a new
  kernel at every bond dimension, but the same number of dispatches;
* padded, jitted -- the composites fuse, so one kernel replaces each of
  those groups. This is the combination Sec. 9 predicted would matter;
* unpadded, jitted -- included precisely because it is the bad case:
  jax.jit traces per input shape, so without padding the compile count
  can rise instead of falling. Printed, not hidden.

Both calculations that the port covers are swept, because they stress
different parts of it: a ground state (contraction chain + Gram/eigh
truncation + Lanczos reorthogonalization) and a real-time TDVP quench
(the Krylov propagator, whose chi grows with time by construction).

Cold and warm are reported separately for the same reason
port_speedup.py does it: compilation is a one-time tax per shape, so a
one-shot script and a sweep inside one process see very different
numbers.

    python3 jit_speedup.py --backend jax --n 8 --maxm 24
    python3 jit_speedup.py --backend jax --cases padded-jit,padded-eager

Pin BLAS threads (MKL_NUM_THREADS=1 OMP_NUM_THREADS=1) before running,
GPU included -- see CLAUDE.md's benchmarking section.
"""

import argparse
import json
import os
import sys
import time


def _add_src(path):
    sys.path.insert(0, path or os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))


def build_chain(n, maxm, nsweeps, seed):
    """A Heisenberg chain. Same builder and same seed for every timing, so
    a difference between two rows is the knob and not the model (DMRG
    starts from a random MPS)."""
    import numpy as np
    from dmrgpy import spinchain
    np.random.seed(seed)
    sc = spinchain.Spin_Chain(["S=1/2" for _ in range(n)],
                              itensor_version="python")
    sc.set_hamiltonian(heisenberg(sc, n))
    sc.maxm = maxm
    sc.nsweeps = nsweeps
    return sc


def heisenberg(sc, n):
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1]
        h = h + sc.Sy[i] * sc.Sy[i + 1]
        h = h + sc.Sz[i] * sc.Sz[i + 1]
    return h


def staggered_field(sc, n):
    h = 0
    for i in range(n):
        h = h + (-1) ** i * sc.Sz[i]
    return h


def run_groundstate(args):
    sc = build_chain(args.n, args.maxm, args.nsweeps, args.seed)
    return float(sc.gs_energy(mode="DMRG").real)


def run_tdvp(args):
    """Real-time quench: ground state of a staggered field, evolved under
    Heisenberg, measuring <Sz_0>(t). Returns the last value so a wrong
    answer cannot pass as a fast one.

    One chain object throughout, with the Hamiltonian swapped after the
    ground state is prepared: an MPS carries its chain's own Index
    identities, so evolving it on a second, independently built chain
    fails outright rather than quietly giving a wrong number."""
    import numpy as np
    from dmrgpy import timedependent
    sc = build_chain(args.n, args.maxm, args.nsweeps, args.seed)
    sc.set_hamiltonian(staggered_field(sc, args.n))
    wf = sc.get_gs()
    sc.set_hamiltonian(heisenberg(sc, args.n))
    _ts, obs = timedependent.evolve_and_measure(sc, operator=sc.Sz[0],
                                                nt=args.nt, dt=args.dt, wf=wf)
    return float(np.real(obs[-1]))


CASES = {
    #  name            pad     jit
    "eager": (False, False),
    "padded-eager": (True, False),
    "padded-jit": (True, True),
    "jit": (False, True),
}


def measure(fn, args, bk, pad, jit, reps):
    """(cold, warm, compilations, value). Cold is the first call in a
    freshly-cleared jit/compile state; warm is the best of `reps` further
    calls, i.e. steady state with every kernel already compiled."""
    bk.set_backend(args.backend)
    bk.set_pad_bonds(args.maxm if pad else None)
    bk.set_jit(jit)
    t0 = time.time()
    value = fn(args)
    cold = time.time() - t0
    warm = None
    for _ in range(reps):
        t0 = time.time()
        value = fn(args)
        dt = time.time() - t0
        warm = dt if warm is None else min(warm, dt)
    comps = bk.compilations()
    bk.set_backend("numpy")
    bk.set_pad_bonds(None)
    bk.set_jit("auto")
    return cold, warm, comps, value


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--src", default=None, help="path to dmrgpy's src/")
    ap.add_argument("--backend", default="jax", choices=("jax", "numpy"),
                    help="array backend to time (numpy ignores both knobs, "
                         "which is itself worth confirming)")
    ap.add_argument("--n", type=int, default=8, help="chain length")
    ap.add_argument("--maxm", type=int, default=24,
                    help="bond dimension, and the padding dimension when "
                         "padding is on")
    ap.add_argument("--nsweeps", type=int, default=3)
    ap.add_argument("--nt", type=int, default=10, help="TDVP time steps")
    ap.add_argument("--dt", type=float, default=0.05)
    ap.add_argument("--reps", type=int, default=2,
                    help="warm repetitions; the best is reported")
    ap.add_argument("--seed", type=int, default=5)
    ap.add_argument("--cases", default=",".join(CASES),
                    help="comma-separated subset of: " + ", ".join(CASES))
    ap.add_argument("--calculations", default="groundstate,tdvp")
    ap.add_argument("--json", default=None, help="write raw results here")
    args = ap.parse_args()

    _add_src(args.src)
    from dmrgpy.pyitensor import backend as bk

    fns = {"groundstate": run_groundstate, "tdvp": run_tdvp}
    out = []
    for calc in args.calculations.split(","):
        print("\n== %s (backend=%s, n=%d, maxm=%d) ==" %
              (calc, args.backend, args.n, args.maxm))
        print("%-14s %10s %10s %8s   %s" %
              ("case", "cold [s]", "warm [s]", "kernels", "value"))
        baseline = None
        for case in args.cases.split(","):
            pad, jit = CASES[case]
            cold, warm, comps, value = measure(fns[calc], args, bk, pad, jit,
                                               args.reps)
            if baseline is None:
                baseline = (cold, warm)
            print("%-14s %10.2f %10.2f %8s   %.12f%s" %
                  (case, cold, warm, "-" if comps is None else comps, value,
                   "" if case == args.cases.split(",")[0] else
                   "   (%.2fx cold, %.2fx warm)" % (baseline[0] / cold,
                                                    baseline[1] / warm)))
            out.append(dict(calculation=calc, case=case, backend=args.backend,
                            n=args.n, maxm=args.maxm, cold=cold, warm=warm,
                            compilations=comps, value=value))
    if args.json:
        with open(args.json, "w") as f:
            json.dump(out, f, indent=1)
        print("\nwrote", args.json)


if __name__ == "__main__":
    main()
