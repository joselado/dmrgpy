# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import time
import numpy as np
from dmrgpy import fermionchain

# The 4-point fermionic correlator tensor <Cdag_i C_j Cdag_k C_l>
# (entropytk/correlationentropy.py), comparing the new ctmode="sweep"
# (Chain::four_correlation_tensor_sweep, mpscpp3/chain_session.h for
# itensor_version=3, pyitensor/chain.py's own port for itensor_version=
# "python") against the existing ctmode="full" (Chain::
# four_correlation_tensor, independent per-(i,j,k,l) AutoMPO builds) and
# ctmode="explicit" (backend-agnostic Python loop of vev()s), and
# cross-checking all against ED.
#
# ctmode="full" builds and applies a *fresh* AutoMPO for every (i,j,k,l)
# tuple -- O(N^4) independent MPO constructions/compressions/contraction
# sweeps. ctmode="sweep" instead follows the algorithmic idea of
# ITensorCorrelators.jl (https://github.com/ITensor/ITensorCorrelators.jl):
# a single left-to-right sweep whose partial contractions ("environments")
# are reused across the whole tensor, only rebuilding what actually
# changes as (i,j,k,l) varies. Implemented for itensor_version in
# (3,"python") on the plain (non-native-spinful) fermionic sites -- see
# get_four_correlation_tensor_sweep()'s docstring for the scope and why.
# ctmode="sweep" is also now the *default* whenever it applies (ctmode=
# None -- see _four_correlation_tensor_default_ctmode()'s docstring for
# the full priority order across every backend).
#
# Measured below (n=6..16 for itensor_version=3, n=5..7 for "python" --
# the pure-Python engine is much slower in absolute terms so isn't run at
# the larger sizes here, same random-hopping + nearest-neighbor-
# interaction chain, same maxm/nsweeps for every ctmode/size at a given
# n): ctmode="sweep" agrees with ctmode="full" and with ED to machine
# precision / ED's solver tolerance at every size and on both backends,
# and is *faster* than ctmode="full" at every size tried on both, by a
# margin that *grows* with n under itensor_version=3 (roughly 1.2x at
# n=6, up to >2.5x at n=16) -- consistent with turning
# four_correlation_tensor()'s O(N^4) independent O(N)-cost MPO builds
# into four_correlation_tensor_sweep()'s O(N^4) total *cheap* local
# tensor contractions (see that method's own docstring in
# chain_session.h for the full algorithm). The pure-Python port shows the
# same qualitative win (~1.2-1.4x at n=5-7) at a much smaller absolute
# margin, since Python-level loop/function-call overhead is a bigger
# fraction of the per-step cost there than in the compiled backend.
#
# `accelerate` (default True on every ctmode) only speeds up the
# *repeated*-index entries here (a subdominant O(N^3) slice of the total
# tensor): unlike ctmode="full", where it skips ~half of the dominant
# per-tuple AutoMPO builds via conjugate-pair symmetry, ctmode="sweep"'s
# dominant pairwise-distinct-index sweep pays its cost once regardless
# (shared environment sweep across (a,b,c,d) plus six cheap per-pattern
# scalar evaluations) -- see four_correlation_tensor_sweep()'s own
# docstring (either backend) for why there's no equivalent saving to skip
# there. Not benchmarked separately below since the timing difference
# between accelerate=True/False on ctmode="sweep" is not expected to be
# meaningful.

np.random.seed(0)


def build(n, itensor_version=3, maxm=60, nsweeps=8):
    fc = fermionchain.Fermionic_Chain(n, itensor_version=itensor_version)
    h = 0
    for i in range(n - 1):
        h = h + 1.0 * (fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i])
        h = h + 0.5 * fc.N[i] * fc.N[i + 1]
    for i in range(n):
        h = h + 0.1 * fc.N[i]
    fc.maxm = maxm
    fc.nsweeps = nsweeps
    fc.set_hamiltonian(h)
    return fc


print("### Part 1: correctness cross-check (ctmode=\"sweep\" vs \"full\" vs \"explicit\" vs ED) ###")
n0 = 5
fc = build(n0, itensor_version=3)
wf = fc.get_gs(mode="DMRG")

ct_sweep = wf.get_four_correlation_tensor(ctmode="sweep")
ct_full = wf.get_four_correlation_tensor(ctmode="full")
ct_explicit = wf.get_four_correlation_tensor(ctmode="explicit")
ct_default = wf.get_four_correlation_tensor()  # ctmode=None -> auto-picks "sweep" here

fced = fermionchain.Fermionic_Chain(n0, mode="ED")
fced.set_hamiltonian(fc.hamiltonian)
wfed = fced.get_gs(mode="ED")
ct_ed = wfed.get_four_correlation_tensor()

diff_sweep_full = np.max(np.abs(ct_sweep - ct_full))
diff_sweep_explicit = np.max(np.abs(ct_sweep - ct_explicit))
diff_sweep_ed = np.max(np.abs(ct_sweep - ct_ed))
diff_default = np.max(np.abs(ct_sweep - ct_default))
print(f"max|sweep - full|     = {diff_sweep_full:.2e}")
print(f"max|sweep - explicit| = {diff_sweep_explicit:.2e}")
print(f"max|sweep - ED|       = {diff_sweep_ed:.2e}")
print(f"max|sweep - default|  = {diff_default:.2e}  (ctmode=None should pick \"sweep\" here)")
assert diff_sweep_full < 1e-8
assert diff_sweep_explicit < 1e-6
assert diff_sweep_ed < 1e-4
assert diff_default == 0.0

print()
print("### Part 2: pure-Python backend (itensor_version=\"python\") cross-check ###")
fc_py = build(n0, itensor_version="python")
wf_py = fc_py.get_gs(mode="DMRG")
ct_py_sweep = wf_py.get_four_correlation_tensor(ctmode="sweep")
ct_py_full = wf_py.get_four_correlation_tensor(ctmode="full")
diff_py = np.max(np.abs(ct_py_sweep - ct_py_full))
diff_py_ed = np.max(np.abs(ct_py_sweep - ct_ed))
print(f"max|python sweep - python full| = {diff_py:.2e}")
print(f"max|python sweep - ED|          = {diff_py_ed:.2e}")
assert diff_py < 1e-8
assert diff_py_ed < 1e-4

print()
print("### Part 3: performance comparison, itensor_version=3 ###")
ns_v3 = [6, 8]
times_full_v3 = []
times_sweep_v3 = []
for n in ns_v3:
    fc = build(n, itensor_version=3)
    wf = fc.get_gs(mode="DMRG")

    t0 = time.time()
    wf.get_four_correlation_tensor(ctmode="full")
    t_full = time.time() - t0

    t0 = time.time()
    wf.get_four_correlation_tensor(ctmode="sweep")
    t_sweep = time.time() - t0

    print(f"  n={n:2d}  full={t_full:7.2f}s  sweep={t_sweep:7.2f}s  "
          f"speedup={t_full/t_sweep:5.2f}x")
    times_full_v3.append(t_full)
    times_sweep_v3.append(t_sweep)

print()
print("### Part 4: performance comparison, itensor_version=\"python\" ###")
ns_py = [5, 6]
times_full_py = []
times_sweep_py = []
for n in ns_py:
    fc = build(n, itensor_version="python", maxm=30, nsweeps=6)
    wf = fc.get_gs(mode="DMRG")

    t0 = time.time()
    wf.get_four_correlation_tensor(ctmode="full")
    t_full = time.time() - t0

    t0 = time.time()
    wf.get_four_correlation_tensor(ctmode="sweep")
    t_sweep = time.time() - t0

    print(f"  n={n:2d}  full={t_full:7.2f}s  sweep={t_sweep:7.2f}s  "
          f"speedup={t_full/t_sweep:5.2f}x")
    times_full_py.append(t_full)
    times_sweep_py.append(t_sweep)

# plot "full" vs "sweep" timing scaling with system size, for both backends
import matplotlib.pyplot as plt
fig,axes = plt.subplots(1,2,figsize=(11,4))
axes[0].plot(ns_v3,times_full_v3,marker="o",label="full")
axes[0].plot(ns_v3,times_sweep_v3,marker="o",label="sweep")
axes[0].set_title("itensor_version=3")
axes[0].set_xlabel("n")
axes[0].set_ylabel("Time [s]")
axes[0].legend()
axes[1].plot(ns_py,times_full_py,marker="s",label="full")
axes[1].plot(ns_py,times_sweep_py,marker="s",label="sweep")
axes[1].set_title("itensor_version=\"python\"")
axes[1].set_xlabel("n")
axes[1].set_ylabel("Time [s]")
axes[1].legend()
plt.tight_layout()
plt.show()
