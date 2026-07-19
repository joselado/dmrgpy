# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Wall-clock timing comparison of all four DMRG dynamical-correlator
# submodes -- KPM, CVM, TD, and the new TDZ (arXiv:2311.10909, see
# tdz.py) -- for a 20-site spin-1/2 Heisenberg chain, computed once on
# the same cached ground state (ITensor v3 backend, verified below --
# if the v3 extension isn't compiled, mode.py silently falls back to ED
# for every call, which is a completely different, non-representative
# cost profile, so this script refuses to mislabel that as "v3"). Unlike
# the v2_VS_v3_* examples in this directory (which cross-check
# *results*), this one is purely about *timing*, like
# backend_timing_gs_energy/.
#
# CVM solves one linear system per requested frequency point (its cost
# is therefore *per point*, not a fixed cost for the whole spectrum like
# the other three submodes). It used to be prohibitively slow to measure
# here (~170.8 s/point, ~1.9 hours for these 40 points, measured before
# cvm.py's CG loop learned to stop at the truncation-imposed residual
# floor instead of always burning the full cvm_nit iteration budget --
# at this system size and cvm_maxm the residual floor is hit almost
# immediately, so nearly all of those iterations changed nothing in the
# returned correction vector, confirmed directly: the pre-fix code's
# value and the post-fix value agree to every printed digit while the
# time per point dropped ~19x). It is now cheap enough to just measure
# live below like the other three submodes. Watch the per-point
# "residual = ..." prints, though: at this size and the default
# cvm_maxm, CG stalls far above cvm_tol, i.e. the CVM numbers time an
# *under-converged* solve -- raising cvm_maxm makes CVM slower but
# actually converged; see docs/user_guide.md's CVM section.
#
# Each get_dynamical_correlator() call below internally re-verifies the
# ground state before computing its own correlator (dynamics.py calls
# set_initial_wf(wf0) unconditionally). Since groundstate.py stopped
# re-sending an unchanged Hamiltonian to the session (which used to
# invalidate the session's energy cache and turn every re-verification
# into a real ~1 s warm re-sweep -- and, for KPM, invalidate the cached
# band edges too), that re-verification is now an actual cache hit, so
# the numbers below are a clean measurement of each submode's own cost.
#
# Confirmed directly on one machine (numbers will vary by machine/load --
# this script is meant to be re-run, not read as a fixed benchmark
# table): at this system size and dt=0.1/delta=0.15 (-> nt =
# int(damping_periods/delta/dt) = 400 internal TDVP time steps, NOT the
# 40 *frequency* points in `es` -- those are just where the single
# resulting time series gets FFT-interpolated onto afterward), TD and
# TDZ have almost identical total wall-clock cost (~49-54s for those 400
# steps). This is expected, not a sign TDZ has no benefit: both are
# capped at the same bond dimension (self.maxm), so once entanglement
# growth pushes either method to that cap, the per-step cost (which
# scales with chi^3) is the same for both -- the paper's actual benefit
# is *accuracy at a given chi* (or equivalently, needing a *smaller* chi
# for a given accuracy) at long times/low frequencies, not raw per-step
# speed at a fixed, already-saturated bond dimension. Seeing that
# benefit requires a fine-delta/long-time accuracy comparison (e.g.
# against KPM or ED), not a pure timing run like this one.
import time
import numpy as np
from dmrgpy import spinchain
from dmrgpy import cppext

n = 20
es = np.linspace(0.0, 3.0, 40)
delta = 0.15

if not cppext.available(3):
    raise RuntimeError(
        "ITensor v3 extension not compiled in this environment -- "
        "every DMRG call below would silently fall back to ED instead "
        "(see mode.py), which is not what this script is meant to time. "
        "Run `python install.py --itensor-version=3` (or `=both`) first.")

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

t0 = time.perf_counter()
sc.get_dynamical_correlator(mode="DMRG", submode="CVM", name=name, es=es, delta=delta)
results["CVM"] = time.perf_counter() - t0

t0 = time.perf_counter()
sc.get_dynamical_correlator(mode="DMRG", submode="TD", name=name, es=es, delta=delta)
results["TD"] = time.perf_counter() - t0

t0 = time.perf_counter()
sc.get_dynamical_correlator(mode="DMRG", submode="TDZ", name=name, es=es, delta=delta,
        alpha0=0.1, n_max=4, dt=0.1)
results["TDZ"] = time.perf_counter() - t0

print(f"{'submode':>8}{'time (s)':>12}   note")
print(f"{'KPM':>8}{results['KPM']:12.2f}   one Chebyshev-moment run, full spectrum")
print(f"{'CVM':>8}{results['CVM']:12.2f}   one linear solve per point x {len(es)} points ({results['CVM']/len(es):.2f} s/point)")
print(f"{'TD':>8}{results['TD']:12.2f}   one real-time evolution run, full spectrum (FFT)")
print(f"{'TDZ':>8}{results['TDZ']:12.2f}   one complex-time evolution run, full spectrum (FFT)")
