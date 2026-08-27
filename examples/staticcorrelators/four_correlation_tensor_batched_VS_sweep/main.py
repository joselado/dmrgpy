# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import time
import numpy as np
from dmrgpy import fermionchain, cppext

# The 4-point fermionic correlator tensor <Cdag_i C_j Cdag_k C_l>, comparing
# ctmode="batched" (pyitensor/fourpoint.py) against the ctmode="sweep" it
# replaces and against the compiled ITensor v3 backend.
#
# What changes. "sweep" and "fold" both evaluate one tuple's worth of MPS
# transfer at a time -- a chain of chi x chi contractions, and there are
# O(n^4) of them. Two things go wrong with that at n=30. The contractions
# are individually far too small to keep BLAS busy, and the *repeated*-index
# tuples (the O(n^3) whose four indices are not pairwise distinct) get no
# environment reuse at all, so they end up costing more than the O(n^4)
# distinct ones: measured 61-65% of the sweep's total runtime at n=12..20.
#
# "batched" reorganizes the same arithmetic around a trie. Written in
# site-sorted order, a tuple becomes a sequence of at most four (local
# matrix, parity) steps, and two tuples that agree on their first r steps
# share a partial environment whatever sites those steps landed on. So all
# n^4 tuples collapse onto a few dozen matrix sequences, each holding one
# array of environments batched over the site combinations that realize it,
# and one site of the sweep is two large GEMMs per node instead of tens of
# thousands of small contractions. Both halves of the tensor -- distinct and
# repeated alike -- go through the same machinery.
#
# The batch axis is also what makes this worth putting on a GPU: see
# docs/gpu_cpu_performance.md, and run with
#     from dmrgpy.pyitensor import backend ; backend.set_backend("jax")
# before the get_four_correlation_tensor call to try it.
#
# Timings printed below are single-threaded -- pin BLAS threads
# (MKL_NUM_THREADS=1 OMP_NUM_THREADS=1) before trusting any of them.

NS = [8, 10, 12, 14, 16]
MAXM = 20


def chain(n, version):
    fc = fermionchain.Fermionic_Chain(n, itensor_version=version)
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i]
        h = h + 1.0 * fc.N[i] * fc.N[i + 1]
    fc.set_hamiltonian(h)
    fc.maxm = MAXM
    fc.nsweeps = 4
    return fc


def timed(wf, **kw):
    t0 = time.time()
    ct = wf.get_four_correlation_tensor(**kw)
    return time.time() - t0, ct


have_v3 = cppext.available(3)
rows = []
for n in NS:
    wf = chain(n, "python").get_gs()
    t_batched, ct_batched = timed(wf, ctmode="batched")
    t_sweep, ct_sweep = timed(wf, ctmode="sweep")
    # the two implementations share no code below the MPS itself, so they
    # agree to machine precision or one of them is wrong
    err = np.abs(ct_batched - ct_sweep).max()
    assert err < 1e-10, "batched disagrees with sweep by %.2e at n=%d" % (err, n)

    t_v3 = np.nan
    if have_v3:
        wf3 = chain(n, 3).get_gs()
        t_v3, _ = timed(wf3)      # v3's own default, i.e. its C++ sweep
    rows.append((n, t_batched, t_sweep, t_v3))
    print("n=%2d  batched %7.2fs   sweep %7.2fs   v3 %7.2fs   agree %.2e"
          % (n, t_batched, t_sweep, t_v3, err), flush=True)

rows = np.array(rows, dtype=float)

import matplotlib.pyplot as plt
fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
ax.plot(rows[:, 0], rows[:, 2], "o-", label='ctmode="sweep" (python)')
if have_v3:
    ax.plot(rows[:, 0], rows[:, 3], "^-", label="ITensor v3 (C++)")
ax.plot(rows[:, 0], rows[:, 1], "s-", label='ctmode="batched" (python)')
ax.set_xlabel("chain length $n$")
ax.set_ylabel("time for the full tensor [s]")
ax.set_yscale("log")
ax.set_title(r"$\langle C^\dagger_i C_j C^\dagger_k C_l\rangle$, maxm=%d" % MAXM)
ax.legend()

ax2.plot(rows[:, 0], rows[:, 2] / rows[:, 1], "s-", label="vs sweep (python)")
if have_v3:
    ax2.plot(rows[:, 0], rows[:, 3] / rows[:, 1], "^-", label="vs ITensor v3 (C++)")
ax2.axhline(1.0, color="k", lw=0.8, ls="--")
ax2.set_xlabel("chain length $n$")
ax2.set_ylabel('speedup of ctmode="batched"')
ax2.set_title("batching wins more as the tuple count grows")
ax2.legend()
plt.tight_layout()
plt.show()
