# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# KPM dynamical correlator on the CPU (NumPy) and on the GPU (JAX), on the
# same chain, with the same seeded random starting state.
#
# The pure-Python backend (itensor_version="python") can hold its tensors
# either in host memory or on a device; pyitensor/backend.py is the switch.
# This script measures both on whatever hardware it is run on, checks that
# they agree, and plots the result.
#
# WHAT TO EXPECT, because a speedup is not guaranteed and the reason is
# physical rather than a matter of installation:
#
#   * the device wins on BOND DIMENSION, not on chain length. Every device
#     operation costs ~0.35 ms however small it is, so many small tensors
#     lose and few large ones win.
#   * measured on an NVIDIA H200 against one Xeon core, KPM is 0.13x at
#     kpmmaxm=40, break-even near 120, 3.3x at 160 and 10.8x at 240 --
#     see docs/gpu_cpu_performance.md for the full table.
#   * so at the small KPMMAXM values this script defaults to (it has to
#     finish on a laptop), a GPU is *expected to be slower*, and on a
#     CPU-only machine "jax" runs on the host and is slower still. Raise
#     KPMMAXM to 160-240 on a real GPU node to see the crossover.
#   * the first device run also pays a one-time XLA compile cost (a kernel
#     per tensor shape), which is why the timings below are taken twice and
#     both are reported: "cold" includes compilation, "warm" does not.
#
# Correctness is checked with the exact zeroth-moment sum rule,
# int S_zz(w) dw = <GS|Sz Sz|GS> = 1/4 for a spin-1/2 site, which holds
# whatever chain, coupling or bond dimension was used -- so a timing point
# can never be a silently wrong spectrum.

import time
import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import spinchain
from dmrgpy.pyitensor import backend

N = 8                          # sites (small so this finishes anywhere)
KPMMAXM = [20, 40]             # bond dimensions; try [160, 240] on a GPU node
MAXM = 30                      # ground-state bond dimension
DELTA = 0.15
ES = np.linspace(-8., 8., 1500)   # wide, so the sum rule sees all the weight
SEED = 1234


def chain(kpmmaxm):
    np.random.seed(SEED)        # same random MPS start on both backends
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(N)],
                              itensor_version="python")
    h = 0
    for i in range(N-1):
        h = h + sc.Sx[i]*sc.Sx[i+1]
        h = h + sc.Sy[i]*sc.Sy[i+1]
        h = h + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = MAXM
    sc.kpmmaxm = kpmmaxm
    return sc


def kpm_once(kpmmaxm):
    """One correlator, timed. Returns (seconds, omegas, spectrum)."""
    sc = chain(kpmmaxm)
    sc.gs_energy(mode="DMRG")                       # not part of the timing
    t0 = time.perf_counter()
    x, y = sc.get_dynamical_correlator(mode="DMRG", submode="KPM",
                                        name=(sc.Sz[0], sc.Sz[0]),
                                        es=ES, delta=DELTA)
    return time.perf_counter()-t0, np.array(x), np.array(y).real


backends = ["numpy"]
jax_on_device = False           # is the JAX backend actually on a GPU?
try:
    backend.set_backend("jax")
    info = backend.device_info()
    jax_on_device = "cuda" in info.lower() or "gpu" in info.lower()
    print("JAX backend available, device:", info)
    if not jax_on_device:
        print("  NOTE: this JAX has no GPU device, so it will run on the host.\n"
              "  That is still a useful check -- it exercises the same code path\n"
              "  a device does -- but it is NOT a CPU-vs-GPU comparison, and\n"
              "  eager JAX on the host is simply slower than NumPy. The labels\n"
              "  below say so.")
    backends.append("jax")
except Exception as exc:
    print("JAX not available (%s), plotting the CPU result only" % exc)
backend.set_backend("numpy")

# Label honestly: "jax" alone would invite reading a host-only run as a GPU
# result, which is the one conclusion this script must not encourage.
JAXLABEL = "jax (GPU)" if jax_on_device else "jax (host, no GPU)"
LABEL = {"numpy": "numpy (CPU)", "jax": JAXLABEL}

results = {}
for name in backends:
    backend.set_backend(name)
    results[name] = {"cold": [], "warm": [], "spectra": []}
    for kpmmaxm in KPMMAXM:
        t_cold, x, y = kpm_once(kpmmaxm)          # includes any compilation
        t_warm, _, _ = kpm_once(kpmmaxm)          # shapes already compiled
        weight = np.trapezoid(y, x)
        results[name]["cold"].append(t_cold)
        results[name]["warm"].append(t_warm)
        results[name]["spectra"].append((x, y))
        print("%-6s kpmmaxm=%-4d cold %7.2fs  warm %7.2fs   "
              "int S dw = %.6f (exact 0.25, error %.1e)"
              % (name, kpmmaxm, t_cold, t_warm, weight, abs(weight-0.25)))
backend.set_backend("numpy")      # leave the process as we found it

if "jax" in results:
    print("\nwarm speedup (numpy time / %s time), >1 means %s wins:"
          % (JAXLABEL, JAXLABEL))
    for i, kpmmaxm in enumerate(KPMMAXM):
        print("  kpmmaxm=%-4d %.2fx" % (
            kpmmaxm, results["numpy"]["warm"][i]/results["jax"]["warm"][i]))

fig, axes = plt.subplots(1, 3 if "jax" in results else 2, figsize=(15, 4.3))

# 1. the spectra themselves, to show the two backends agree
ax = axes[0]
for name in backends:
    for kpmmaxm, (x, y) in zip(KPMMAXM, results[name]["spectra"]):
        ax.plot(x, y, ls="-" if name == "numpy" else "--",
                label="%s, kpmmaxm=%d" % (LABEL[name], kpmmaxm))
ax.set_xlim(-1, 6)
ax.set_xlabel("Energy $\\omega$")
ax.set_ylabel("$S_{zz}(\\omega)$")
ax.set_title("same spectrum on both backends")
ax.legend(fontsize=7)

# 2. absolute time, cold and warm
ax = axes[1]
for name in backends:
    ax.plot(KPMMAXM, results[name]["warm"], "o-",
            label="%s warm" % LABEL[name])
    ax.plot(KPMMAXM, results[name]["cold"], "s--", alpha=.6,
            label="%s cold" % LABEL[name])
ax.set_xlabel("kpmmaxm (KPM bond dimension)")
ax.set_ylabel("time for one correlator [s]")
ax.set_yscale("log")
ax.set_title("cost (cold includes XLA compilation)")
ax.legend(fontsize=7)
ax.grid(alpha=.3)

# 3. the speedup, with break-even marked
if "jax" in results:
    ax = axes[2]
    sp = [results["numpy"]["warm"][i]/results["jax"]["warm"][i]
          for i in range(len(KPMMAXM))]
    ax.plot(KPMMAXM, sp, "o-")
    ax.axhline(1.0, color="k", ls="--", lw=1)
    ax.set_xlabel("kpmmaxm (KPM bond dimension)")
    ax.set_ylabel("numpy time / %s time" % JAXLABEL)
    ax.set_title(("above the dashed line the device wins\n"
                  "(crossover is near kpmmaxm~120 on an H200)")
                 if jax_on_device else
                 ("no GPU here: eager JAX on the host is slower than NumPy\n"
                  "by construction -- see docs/gpu_cpu_performance.md"))
    ax.grid(alpha=.3)

fig.tight_layout()
plt.show()
