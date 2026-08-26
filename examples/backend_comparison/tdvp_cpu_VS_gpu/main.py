# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Real-time TDVP evolution on the CPU (NumPy) and on the GPU (JAX), on the
# same chain, from the same seeded ground state -- the time-evolution
# counterpart of examples/backend_comparison/kpm_cpu_VS_gpu.
#
# WHY TIME EVOLUTION IS THE INTERESTING CASE FOR A DEVICE. The pure-Python
# backend's GPU port wins on BOND DIMENSION and on nothing else: every
# device operation costs ~0.35 ms however small it is, so the crossover
# sits at chi ~ 120-160 (see docs/gpu_cpu_performance.md). A 1D ground
# state cannot get there -- an area law with log corrections caps its
# entanglement, and a uniform Heisenberg chain converges at chi ~ 60
# whatever maxm says. Real-time evolution has no such cap: entanglement
# grows roughly linearly in time after a quench, so chi climbs into the
# paying range by construction, simply by evolving for longer.
#
# WHAT TO EXPECT HERE. The defaults are small enough to finish in a couple
# of minutes on a laptop -- eager JAX on a host is many times slower than
# NumPy at this size, so they are chosen for the *slowest* configuration
# below rather than the fastest. They therefore sit far *below* the
# crossover: on a real GPU node this
# script is expected to be slower than NumPy until NT and MAXM are raised
# (try MAXM=240, NT=200), and on a CPU-only machine "jax" runs on the host
# and is slower still. The labels below say which of the two you got.
#
# The third configuration, "jax + jit + pad", is the dispatch-floor knob
# pair: backend.set_pad_bonds freezes every tensor shape so XLA stops
# compiling a new kernel per bond dimension, and backend.set_jit fuses the
# engine's hot inner kernels (contraction, Gram+eigh truncation, the
# Krylov recurrence) into one compiled kernel each. Neither changes what is
# computed -- which is exactly what the middle panel checks.
#
# Correctness is anchored on exact diagonalization: the same quench is
# evolved with mode="ED", so a fast trajectory can never be a silently
# wrong one.

import time
import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import spinchain
from dmrgpy import timedependent
from dmrgpy.pyitensor import backend

N = 6                      # sites (small so ED is possible and this finishes)
MAXM = 16                  # bond dimension; raise to 240 on a GPU node
NT = 12                    # time steps; raise to 200 on a GPU node
DT = 0.05
SEED = 1234


def chain():
    """The chain, with the same seed every time: DMRG starts from a random
    MPS, so a cross-backend comparison is only meaningful from the same
    start."""
    np.random.seed(SEED)
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(N)],
                              itensor_version="python")
    sc.maxm = MAXM
    sc.nsweeps = 4        # enough for the initial state; the prep is untimed
    return sc


def staggered_field(sc):
    h = 0
    for i in range(N):
        h = h + (-1)**i*sc.Sz[i]          # Neel-favoring field
    return h


def heisenberg(sc):
    h = 0
    for i in range(N-1):
        h = h + sc.Sx[i]*sc.Sx[i+1]
        h = h + sc.Sy[i]*sc.Sy[i+1]
        h = h + sc.Sz[i]*sc.Sz[i+1]
    return h


def prepare(mode="DMRG"):
    """The quench's initial state: the ground state of the staggered
    field, on a chain whose Hamiltonian is then swapped for Heisenberg.

    One chain object throughout: an MPS carries its own chain's Index
    identities, so evolving it on a second, independently built chain
    fails outright. Untimed on purpose -- this is DMRG, and what this
    script is about is the evolution."""
    sc = chain()
    sc.set_hamiltonian(staggered_field(sc))
    wf = sc.get_gs(mode=mode)
    sc.set_hamiltonian(heisenberg(sc))
    return sc, wf


def evolve(sc, wf, mode="DMRG"):
    """One timed trajectory. Returns (seconds, ts, sz)."""
    t0 = time.perf_counter()
    ts, sz = timedependent.evolve_and_measure(sc, operator=sc.Sz[0],
                                              nt=NT, dt=DT, wf=wf, mode=mode)
    return time.perf_counter()-t0, np.array(ts).real, np.array(sz).real


# --- reference: exact diagonalization -----------------------------------
_sc, _wf = prepare(mode="ED")
_t, ts_ed, sz_ed = evolve(_sc, _wf, mode="ED")
print("ED reference done (%d time steps)" % NT)

# --- the configurations to compare --------------------------------------
jax_on_device = False
configs = [("numpy (CPU)", "numpy", None, "auto")]
try:
    backend.set_backend("jax")
    info = backend.device_info()
    jax_on_device = "cuda" in info.lower() or "gpu" in info.lower()
    print("JAX backend available, device:", info)
    if not jax_on_device:
        print("  NOTE: this JAX has no GPU device, so it runs on the host.\n"
              "  That still exercises the same code path a device does --\n"
              "  same immutability, same dispatch, same host round trips --\n"
              "  but it is NOT a CPU-vs-GPU comparison, and eager JAX on the\n"
              "  host is simply slower than NumPy. The labels say so.")
    tag = "GPU" if jax_on_device else "host, no GPU"
    configs.append(("jax eager (%s)" % tag, "jax", None, False))
    configs.append(("jax + jit + pad (%s)" % tag, "jax", MAXM, True))
except Exception as exc:
    print("JAX not available (%s), plotting the CPU result only" % exc)
backend.set_backend("numpy")

results = []
for label, name, pad, jit in configs:
    backend.set_backend(name)
    backend.set_pad_bonds(pad)
    backend.set_jit(jit)
    sc, wf = prepare()                      # untimed: DMRG, not the point
    t_cold, ts, sz = evolve(sc, wf)         # includes any XLA compilation
    t_warm, _, sz2 = evolve(sc, wf)         # shapes already compiled
    assert np.max(np.abs(sz2-sz)) < 1e-12   # evolving twice from the same
    err = np.max(np.abs(sz-sz_ed))          # state must give one answer
    results.append(dict(label=label, ts=ts, sz=sz, cold=t_cold, warm=t_warm,
                        err=err, kernels=backend.compilations()))
    print("%-28s cold %7.2fs  warm %7.2fs   max|<Sz>-ED| = %.2e   kernels=%s"
          % (label, t_cold, t_warm, err, backend.compilations()))
backend.set_backend("numpy")     # leave the process as we found it
backend.set_pad_bonds(None)
backend.set_jit("auto")

# every configuration must reproduce the exact trajectory
for r in results:
    assert r["err"] < 1e-6, "%s disagrees with ED by %.2e" % (r["label"], r["err"])
print("\nall configurations agree with exact diagonalization")

# --- plot ----------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))

ax = axes[0]
ax.plot(ts_ed, sz_ed, "k-", lw=3, alpha=.35, label="exact diagonalization")
for r in results:
    ax.plot(r["ts"], r["sz"], "--", lw=1.4, label=r["label"])
ax.set_xlabel("time")
ax.set_ylabel(r"$\langle S_z^0\rangle(t)$")
ax.set_title("TDVP quench, n=%d, maxm=%d" % (N, MAXM))
ax.legend(fontsize=8)

ax = axes[1]
for r in results[1:]:
    ax.semilogy(r["ts"], np.abs(r["sz"]-results[0]["sz"])+1e-18, label=r["label"])
ax.axhline(1e-10, color="k", ls=":", lw=1)
ax.set_xlabel("time")
ax.set_ylabel(r"$|\langle S_z^0\rangle - \langle S_z^0\rangle_{\rm numpy}|$")
ax.set_title("device vs host: the same trajectory\n(the knobs are pure performance)")
if len(results) > 1:
    ax.legend(fontsize=8)

ax = axes[2]
idx = np.arange(len(results))
ax.bar(idx-0.2, [r["cold"] for r in results], 0.4, label="cold (with compile)")
ax.bar(idx+0.2, [r["warm"] for r in results], 0.4, label="warm")
ax.set_xticks(idx)
ax.set_xticklabels([r["label"].replace(" (", "\n(") for r in results], fontsize=7)
ax.set_ylabel("time for %d TDVP steps [s]" % NT)
ax.set_title("cost (this machine; the crossover\nis at chi~120-160, see the header)")
ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig("tdvp_cpu_VS_gpu.png", dpi=120)
plt.show()
