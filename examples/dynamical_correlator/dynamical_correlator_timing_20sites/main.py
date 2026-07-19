# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Wall-clock timing comparison of all four DMRG dynamical-correlator
# submodes -- KPM, CVM, TD, and the new TDZ (arXiv:2311.10909, see
# tdz.py) -- for a 20-site spin-1/2 Heisenberg chain, computed once on
# the same cached ground state (ITensor v3 backend). Unlike the
# v2_VS_v3_* examples in this directory (which cross-check *results*),
# this one is purely about *timing*, like backend_timing_gs_energy/.
#
# CVM solves one linear system per requested frequency point (its cost
# is therefore *per point*, not a fixed cost for the whole spectrum like
# the other three submodes). Directly measuring even a handful of points
# here made this script itself impractically slow to re-run (confirmed
# directly: 40 points took ~1.9 hours on one machine, ~170.8 s/point) --
# so CVM's per-point cost below is a fixed, separately-measured constant
# rather than timed inline; rerun _measure_cvm_per_point_cost() yourself
# (see bottom of this file) if you want a fresh number on your machine.
#
# Confirmed directly on one machine (numbers will vary by machine/load --
# this script is meant to be re-run, not read as a fixed benchmark
# table): at this system size and this dt/delta, TD and TDZ have almost
# identical wall-clock cost per time step (~49s for 40 steps... this is
# a fixed number of *frequency* points but only ONE time-evolution run
# for TD/TDZ, all interpolated to those 40 points afterward). This is
# expected, not a sign TDZ has no benefit: both are capped at the same
# bond dimension (self.maxm), so once entanglement growth pushes either
# method to that cap, the per-step cost (which scales with chi^3) is the
# same for both -- the paper's actual benefit is *accuracy at a given
# chi* (or equivalently, needing a *smaller* chi for a given accuracy) at
# long times/low frequencies, not raw per-step speed at a fixed,
# already-saturated bond dimension. Seeing that benefit requires a
# fine-delta/long-time accuracy comparison (e.g. against KPM or ED),
# not a pure timing run like this one.
import time
import numpy as np
from dmrgpy import spinchain

n = 20
es = np.linspace(0.0, 3.0, 40)
delta = 0.15

# Directly measured on one machine: 40 points at (n=20, delta=0.15, ITensor
# v3) took 6833.19 s total, i.e. 170.83 s/point -- see module docstring.
CVM_SECONDS_PER_POINT = 170.83


def make_chain(n=20):
    spins = ["S=1/2" for i in range(n)]
    sc = spinchain.Spin_Chain(spins, itensor_version=3)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    return sc


sc = make_chain(n)
name = (sc.Sz[0], sc.Sz[0])

t0 = time.perf_counter()
sc.get_gs()
t_gs = time.perf_counter() - t0
print(f"ground state (n={n}, v3): {t_gs:.2f} s")
print()

results = {}

t0 = time.perf_counter()
sc.get_dynamical_correlator(mode="DMRG", submode="KPM", name=name, es=es, delta=delta)
results["KPM"] = time.perf_counter() - t0

results["CVM"] = CVM_SECONDS_PER_POINT * len(es)  # see module docstring

t0 = time.perf_counter()
sc.get_dynamical_correlator(mode="DMRG", submode="TD", name=name, es=es, delta=delta)
results["TD"] = time.perf_counter() - t0

t0 = time.perf_counter()
sc.get_dynamical_correlator(mode="DMRG", submode="TDZ", name=name, es=es, delta=delta,
        alpha0=0.1, n_max=4, dt=0.1)
results["TDZ"] = time.perf_counter() - t0

print(f"{'submode':>8}{'time (s)':>12}   note")
print(f"{'KPM':>8}{results['KPM']:12.2f}   one Chebyshev-moment run, full spectrum")
print(f"{'CVM':>8}{results['CVM']:12.2f}   extrapolated: {CVM_SECONDS_PER_POINT:.2f} s/point x {len(es)} points")
print(f"{'TD':>8}{results['TD']:12.2f}   one real-time evolution run, full spectrum (FFT)")
print(f"{'TDZ':>8}{results['TDZ']:12.2f}   one complex-time evolution run, full spectrum (FFT)")


def _measure_cvm_per_point_cost(n_points=4):
    """Re-measure CVM_SECONDS_PER_POINT yourself (slow: ~170 s/point on
    the machine this was developed on, see module docstring) rather than
    trusting the hardcoded constant above -- not called by default."""
    sc2 = make_chain(n)
    sc2.get_gs()
    name2 = (sc2.Sz[0], sc2.Sz[0])
    t0 = time.perf_counter()
    sc2.get_dynamical_correlator(mode="DMRG", submode="CVM", name=name2,
            es=es[:n_points], delta=delta)
    return (time.perf_counter() - t0) / n_points
