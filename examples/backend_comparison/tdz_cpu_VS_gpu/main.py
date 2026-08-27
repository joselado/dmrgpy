# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# The complex-time-evolution dynamical correlator (submode="TDZ",
# arXiv:2311.10909) on the CPU (NumPy) and on the GPU (JAX), and the price
# of each of the two optimizations that make it device-viable -- the
# complex-time counterpart of examples/backend_comparison/tdvp_cpu_VS_gpu.
#
# WHAT MAKES TDZ DIFFERENT FROM THE OTHER GPU EXAMPLES. Everywhere else in
# this directory the device story is bond dimension: an operation costs
# ~0.35 ms on a device however small it is, so the crossover sits at
# chi ~ 120-160 (docs/gpu_cpu_performance.md). TDZ is the calculation that
# deliberately works in the *opposite* regime -- the whole point of the
# complex-time contour is that it damps high-energy content, so the bond
# dimension needed for a given accuracy is small (the paper reports
# chi ~ 20-30, against 500-700 for real-time evolution). It therefore sits
# on the wrong side of that crossover by construction, and this script is
# honest about that: expect NumPy to win at the defaults below.
#
# WHAT WAS OPTIMIZED, AND WHY IT IS NOT ABOUT FLOPS. On a device the cost
# of a long chain of small operations is dominated by *synchronization*,
# not arithmetic. JAX dispatch is asynchronous -- the host runs ahead
# enqueueing kernels while the device works -- and every value read back to
# the host drains that queue. Two changes attack exactly that, and this
# script prices them separately:
#
#   krylov  the Krylov exponentiator inside TDVP used to read two numbers
#           (the Lanczos alpha and beta) back to the host per iteration:
#           ~1000 pipeline stalls per time step at n=30. It now keeps them
#           on the device and brings a block home per checkpoint, running
#           speculatively past its own stopping point and then rolling back
#           to the exact same Krylov dimension -- so the answer does not
#           change, which is what the middle panel checks.
#   bras    the n_max+1 phi^(n) overlaps each step needs are taken as one
#           batched sweep instead of n_max+1 separate ones. The bras are
#           fixed for the whole run, so they are stacked and padded once.
#
# Neither changes what is computed. That is the claim the middle panel
# exists to falsify, and it is checked with a real assert as well as
# plotted, so this doubles as a regression test in the spirit of
# examples/dynamical_correlator/dynamical_correlator_tdz_VS_ED.

import time
import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import spinchain
from dmrgpy import tdz as tdzmod
from dmrgpy.pyitensor import backend
from dmrgpy.pyitensor import tdvp

N = 8                      # sites (small so this finishes on a laptop)
MAXM = 16                  # bond dimension; raise to 240 on a GPU node
NT = 20                    # complex-time steps; raise to 200 on a GPU node
DT = 0.1
SEED = 1234
ES = np.linspace(-0.5, 4.0, 200)


def chain():
    """The chain, with the same seed every time: DMRG starts from a random
    MPS, so a cross-backend comparison is only meaningful from the same
    start."""
    np.random.seed(SEED)
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(N)],
                              itensor_version="python")
    h = 0
    for i in range(N-1):
        h = h + sc.Sx[i]*sc.Sx[i+1]
        h = h + sc.Sy[i]*sc.Sy[i+1]
        h = h + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = MAXM
    sc.nsweeps = 3
    return sc


def prepare():
    """The chain with its ground state already solved, untimed.

    Solving it inside the timed region would be defensible -- it is part
    of what one correlator costs -- but it is the same cost in every
    configuration here and, on eager JAX, comparable to the evolution
    itself, so it would dilute exactly the difference this script is
    about. benchmarks/gpu/tdz_bench.py makes the other choice and reports
    the two separately."""
    sc = chain()
    sc.get_gs()               # cached, so the timed call below reuses it
    return sc


def tdz_run(sc):
    """One timed TDZ correlator on an already-prepared chain."""
    t0 = time.perf_counter()
    es, spec = sc.get_dynamical_correlator(
            mode="DMRG", submode="TDZ", name=(sc.Sz[N//2], sc.Sz[N//2]),
            dt=DT, nt=NT, es=ES)
    return time.perf_counter()-t0, np.array(es).real, np.array(spec)


# --- the configurations to compare --------------------------------------
# (label, backend, pad, jit, defer_krylov_sync, batched_bras)
jax_on_device = False
configs = [("numpy, both opts", "numpy", None, "auto", None, True)]
try:
    backend.set_backend("jax")
    info = backend.device_info()
    jax_on_device = "cuda" in info.lower() or "gpu" in info.lower()
    print("JAX backend available, device:", info)
    if not jax_on_device:
        print("  NOTE: this JAX has no GPU device, so it runs on the host.\n"
              "  That exercises the same code path a device does -- same\n"
              "  immutability, same dispatch, same host round trips -- but\n"
              "  the synchronizations this script is about are nearly free\n"
              "  on a host, so the 'krylov' row will look like noise here.\n"
              "  It is NOT a CPU-vs-GPU comparison. The labels say so.")
    tag = "GPU" if jax_on_device else "host, no GPU"
    # Forced off/on rather than left to the automatic device default, so
    # each row prices one change against the same baseline.
    configs += [
        ("jax, neither (%s)" % tag,   "jax", MAXM, True, False, False),
        ("jax, +krylov (%s)" % tag,   "jax", MAXM, True, True,  False),
        ("jax, +bras (%s)" % tag,     "jax", MAXM, True, False, True),
        ("jax, both (%s)" % tag,      "jax", MAXM, True, True,  True),
    ]
except Exception as exc:
    print("JAX not available (%s), plotting the CPU result only" % exc)
backend.set_backend("numpy")

results = []
for label, name, pad, jit, defer, batched in configs:
    backend.set_backend(name)
    backend.set_pad_bonds(pad)
    backend.set_jit(jit)
    tdvp.set_krylov_defer_sync(defer)
    tdzmod.set_batched_bras(batched)
    sc = prepare()                      # untimed: DMRG, not the point
    t_cold, es, spec = tdz_run(sc)      # includes any XLA compilation
    t_warm, _, spec2 = tdz_run(sc)      # shapes already compiled
    assert np.max(np.abs(spec2-spec)) < 1e-12   # same input, same answer
    results.append(dict(label=label, es=es, spec=spec, cold=t_cold,
                        warm=t_warm, kernels=backend.compilations()))
    print("  %-28s cold %7.2fs  warm %7.2fs  kernels %s"
          % (label, t_cold, t_warm, backend.compilations()))
    sys.stdout.flush()

backend.set_backend("numpy")
backend.set_pad_bonds(None)
backend.set_jit("auto")
tdvp.set_krylov_defer_sync(None)
tdzmod.set_batched_bras(True)

# --- the claim: none of this changed the physics ------------------------
ref = results[0]["spec"]
devs = [float(np.max(np.abs(r["spec"]-ref))) for r in results]
for r, d in zip(results, devs):
    assert d < 1e-8, "%s changed the spectrum by %.2e" % (r["label"], d)
print("\nevery configuration agrees with NumPy to %.2e" % max(devs))

# --- plot ----------------------------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))

for r in results:
    ax[0].plot(r["es"], r["spec"].real, lw=1.6, alpha=0.75, label=r["label"])
ax[0].set_xlabel("$\\omega$")
ax[0].set_ylabel("$S_{zz}(\\omega)$")
ax[0].set_title("TDZ correlator: every configuration overlaid")
ax[0].legend(fontsize=7)

ax[1].semilogy(range(len(results)), [max(d, 1e-18) for d in devs], "o-")
ax[1].set_xticks(range(len(results)))
ax[1].set_xticklabels([r["label"] for r in results], rotation=30,
                      ha="right", fontsize=7)
ax[1].axhline(1e-8, color="r", ls="--", lw=1, label="assert threshold")
ax[1].set_ylabel("max deviation from NumPy")
ax[1].set_title("the optimizations change no number")
ax[1].legend(fontsize=7)

x = np.arange(len(results))
ax[2].bar(x-0.2, [r["cold"] for r in results], 0.4, label="cold (1st run)")
ax[2].bar(x+0.2, [r["warm"] for r in results], 0.4, label="warm (2nd run)")
ax[2].set_xticks(x)
ax[2].set_xticklabels([r["label"] for r in results], rotation=30,
                      ha="right", fontsize=7)
ax[2].set_ylabel("seconds")
ax[2].set_title("wall time (N=%d, maxm=%d, nt=%d)" % (N, MAXM, NT))
ax[2].legend(fontsize=7)

plt.tight_layout()
plt.show()
