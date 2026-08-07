# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Compare wall-clock time for the dynamical METTS finite-temperature
# correlator (Many_Body_Chain.metts_dynamical_correlator, see
# vevtk/mettsdynamicalcorrelator.py and pyitensor/metts.py) between ITensor
# v3 (mpscpp3) and the pure-Python backend (itensor_version="python").
# Unlike examples/groundstate/backend_timing_gs_energy (a handful of DMRG
# sweeps), this workload calls the same TDVP effective-Hamiltonian matvec
# thousands of times more per run (many real-time steps x many Lanczos/
# Krylov iterations x many METTS samples) -- the regime where pyitensor
# used to be dramatically slower than v3, not just modestly slower, because
# make_matvec()'s old default fallback re-derived its whole contraction
# plan from scratch on every single matvec call (cheap enough per call to
# be invisible for a DMRG sweep, but the dominant cost here). See
# pyitensor/kernels.py's module docstring ("plain path" section) for the
# fix and the profiling that found it.
#
# Confirmed directly on one machine (numbers will vary by machine/load --
# this script is meant to be re-run, not read as a fixed benchmark table):
# n=10, nsamples=40, nt=10 went from 152.0s (before the kernels.py/tdvp.py
# fixes in this same change) to ~51s afterward -- v3 itself takes ~9s, so
# the remaining gap (~5.7x) is the dense-tensor-vs-QN-block-sparse gap
# real ITensor v3 still has here (unlike ground-state DMRG, TDVP evolves a
# 2-site tensor through many real-time steps, so mpscpp3's ConserveQNs off
# + random-MPS-start tradeoff -- see documentation.md -- costs it less
# than it does for reaching a converged ground state from scratch).
import time
import matplotlib.pyplot as plt
from dmrgpy import spinchain, cppext


def metts_dynamical_correlator_timed(itensor_version, n, nsamples, nt):
    spins = ["S=1/2" for i in range(n)]
    sc = spinchain.Spin_Chain(spins, itensor_version=itensor_version)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    name = (sc.Sz[0], sc.Sz[0])
    t0 = time.perf_counter()
    sc.metts_dynamical_correlator(name, T=1.0, nt=nt, dt=0.1, nsamples=nsamples,
                                   nwarmup=10, dbeta_half_step=0.08, seed=42)
    return time.perf_counter() - t0


# Warm up the python backend once, outside the timed loop, so one-time
# import/first-call costs don't distort the smallest run.
metts_dynamical_correlator_timed("python", 4, nsamples=1, nt=2)

n, nt = 10, 10
sizes_nsamples = [10, 40]

have_v3 = cppext.available(3)
header = f"{'nsamples':>10}"
if have_v3:
    header += f"{'v3 (s)':>10}"
header += f"{'python (s)':>12}"
if have_v3:
    header += f"{'python/v3':>12}"
print(header)

times3, timespy = [], []
for nsamples in sizes_nsamples:
    tpy = metts_dynamical_correlator_timed("python", n, nsamples, nt)
    row = f"{nsamples:10d}"
    if have_v3:
        t3 = metts_dynamical_correlator_timed(3, n, nsamples, nt)
        row += f"{t3:10.2f}"
        times3.append(t3)
    row += f"{tpy:12.2f}"
    timespy.append(tpy)
    if have_v3:
        row += f"{tpy/t3:12.2f}"
    print(row)

if not have_v3:
    print("\n(ITensor v3 not compiled in this environment -- only the "
          "python backend's own timing is shown; run `python install.py "
          "--itensor-version=3` to compile it and see the v3 comparison.)")

plt.figure(figsize=(6, 4.5))
plt.plot(sizes_nsamples, timespy, "o-", label="python")
if have_v3:
    plt.plot(sizes_nsamples, times3, "s-", label="v3 (mpscpp3)")
plt.xlabel("nsamples")
plt.ylabel("wall time [s]")
plt.yscale("log")
plt.title(f"METTS dynamical correlator timing (n={n}, nt={nt})")
plt.legend()
plt.grid(alpha=0.3, which="both")
plt.tight_layout()
plt.show()
