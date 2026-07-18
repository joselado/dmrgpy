# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

# Compare wall-clock time for ground state energy calculations across the
# three DMRG backends -- ITensor v2 (mpscpp2), ITensor v3 (mpscpp3), and
# the pure-Python backend (itensor_version="python", src/dmrgpy/pyitensor/,
# no compiler/pybind11 needed) -- for a spin-1/2 Heisenberg chain, at a few
# system sizes. Unlike the v2_VS_v3_* examples in this directory (which
# cross-check *results* between backends), this one is purely about
# *timing*: how much slower is the compiler-free pure-Python backend
# compared to compiled ITensor, and how does that gap scale with system
# size. If v2/v3 aren't compiled in this environment, their gs_energy()
# calls silently fall back to ED (see mode.py) -- so their timings would
# then reflect ED, not real DMRG; itensor_version="python" always exercises
# the real in-process pyitensor DMRG sweep regardless.
#
# Confirmed directly on one machine (numbers will vary by machine/load --
# this script is meant to be re-run, not read as a fixed benchmark table):
# with numba installed (pyitensor's default accelerated kernel path, see
# pyitensor/kernels.py), the python backend was only ~1.5-2.5x slower than
# v3 for n=8-20, but the gap widened with n (~4-5x by n=28-32) -- v3
# benefits from real ITensor's QN block-sparsity, which pyitensor's dense
# NumPy tensors don't have even with numba JIT acceleration. Without numba,
# expect a much larger gap (see pyitensor/__init__.py's own docstring).
import time
from dmrgpy import spinchain


def gs_energy_timed(itensor_version, n):
    spins = ["S=1/2" for i in range(n)]
    sc = spinchain.Spin_Chain(spins, itensor_version=itensor_version)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    t0 = time.perf_counter()
    energy = sc.gs_energy()
    dt = time.perf_counter() - t0
    return energy, dt


# Warm up the python backend's numba JIT once, outside the timed loop --
# the first call to each numba-compiled kernel pays a one-time compile
# cost (~0.3-0.5s) that has nothing to do with the DMRG algorithm itself
# and would otherwise distort the smallest system size the most.
gs_energy_timed("python", 6)

sizes = [8, 12, 16, 20]
times2, times3, timespy = [], [], []

print(f"{'n':>4}{'v2 (s)':>10}{'v3 (s)':>10}{'python (s)':>12}{'python/v3':>12}")
for n in sizes:
    _, t2 = gs_energy_timed(2, n)
    _, t3 = gs_energy_timed(3, n)
    _, tpy = gs_energy_timed("python", n)
    times2.append(t2) ; times3.append(t3) ; timespy.append(tpy)
    print(f"{n:4d}{t2:10.3f}{t3:10.3f}{tpy:12.3f}{tpy/t3:12.1f}")

print()
print("Total time across all sizes above:")
print(f"  ITensor v2  : {sum(times2):.3f} s")
print(f"  ITensor v3  : {sum(times3):.3f} s")
print(f"  pure Python : {sum(timespy):.3f} s")
