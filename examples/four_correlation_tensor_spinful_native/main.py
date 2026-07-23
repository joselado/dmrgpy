# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import time
import numpy as np
from dmrgpy import fermionchain

# The 4-point fermionic correlator tensor <Cdag_i C_j Cdag_k C_l>
# (entropytk/correlationentropy.py, examples/four_correlation_tensor is the
# spinless version of this script) on the spinful Hubbard chain, comparing
# fermionchain.Spinful_Fermionic_Chain_Native (native Electron/Hubbard sites,
# itensor_version=3 only) against the existing interleaved
# fermionchain.Spinful_Fermionic_Chain. Both classes expose the same flat
# self.C/self.Cdag convention (mode 2*i=up, 2*i+1=down), so the resulting
# (2n,2n,2n,2n) tensors are directly, index-for-index comparable.
#
# Only ctmode="explicit" (a Python loop over vev()s of MultiOperator
# products, backend-agnostic) works for the native class: ctmode="full"
# (mpscpp3/chain_session.h's Chain::four_correlation_tensor, ITensor's own
# AutoMPO + automatic Jordan-Wigner) is hardcoded to the literal "Cdag"/"C"
# operator names, which ITensor's ElectronSite doesn't define (only
# Cup/Cdn/Cdagup/Cdagdn are) -- so the comparison below is native's only
# option against both of Spinful_Fermionic_Chain's options.
#
# Unlike the ground-state/TDVP/KPM benchmarks in
# examples/hubbard_chain_native_sites (where native sites lose to the
# interleaved chain), this one flips: the 4-point tensor is a *static*
# single-shot MPO-MPS-MPS overlap per (i,j,k,l), not an iterative two-site
# variational search, so native sites' "half as many sites, same total
# mode count" advantage shows through without the offsetting two-site
# combined-local-dimension penalty that dominates ground-state DMRG.
#
# That win is size-bounded, not asymptotic, though: Part 2 below only
# sweeps n=3..6 (fast, growing margin in native's favor there), but a
# separate n=12 check (24 flat modes; not included in the loop below
# since it alone takes ~12 minutes per class) found the two back to
# essentially tied (~700s each, native's explicit slightly *behind*
# doubled's C++-accelerated full at 707s vs 697s) -- the O(n^4) growth
# in the number of (i,j,k,l) overlaps evaluated eventually swamps
# native's per-overlap advantage.

n = 4
U = 1.3


def build(chain_cls, maxm=40, nsweeps=14):
    fc = chain_cls(n)
    h = 0
    for i in range(n - 1):
        h = h + (0.9 + 0.2j) * fc.Cdagup[i] * fc.Cup[i + 1]
        h = h + (0.9 + 0.2j) * fc.Cdagdn[i] * fc.Cdn[i + 1]
    for i in range(n):
        h = h + U * (fc.Nup[i] - .5) * (fc.Ndn[i] - .5)
    h = h + h.get_dagger()
    h = h + 0.3 * (fc.Nup[0] - fc.Ndn[0])  # lifts the up/down GS degeneracy
    fc.maxm = maxm
    fc.nsweeps = nsweeps
    fc.set_hamiltonian(h)
    fc.setup_cpp(3)
    return fc


print("### Part 1: correctness cross-check ###")
fc_native = build(fermionchain.Spinful_Fermionic_Chain_Native)
wf_native = fc_native.get_gs(mode="DMRG")
ct_native = wf_native.get_four_correlation_tensor(ctmode="explicit")

fc_double = build(fermionchain.Spinful_Fermionic_Chain)
wf_double = fc_double.get_gs(mode="DMRG")
ct_double = wf_double.get_four_correlation_tensor(ctmode="full", accelerate=True)

diff = np.max(np.abs(ct_native - ct_double))
print(f"max|native - doubled| = {diff:.2e}")
assert diff < 1e-4

fc_ed = build(fermionchain.Spinful_Fermionic_Chain_Native)
wf_ed = fc_ed.get_gs(mode="ED")
ct_ed = wf_ed.get_four_correlation_tensor()
diff_ed = np.max(np.abs(ct_native - ct_ed))
print(f"max|native - ED|      = {diff_ed:.2e}")
assert diff_ed < 1e-4

print()
print("### Part 2: performance comparison ###")
for n_size in [3, 4, 5, 6]:
    n = n_size
    fc_native = build(fermionchain.Spinful_Fermionic_Chain_Native)
    wf_native = fc_native.get_gs(mode="DMRG")
    t0 = time.time()
    wf_native.get_four_correlation_tensor(ctmode="explicit")
    t_native = time.time() - t0

    fc_double = build(fermionchain.Spinful_Fermionic_Chain)
    wf_double = fc_double.get_gs(mode="DMRG")
    t0 = time.time()
    wf_double.get_four_correlation_tensor(ctmode="explicit")
    t_double_explicit = time.time() - t0
    t0 = time.time()
    wf_double.get_four_correlation_tensor(ctmode="full", accelerate=True)
    t_double_full = time.time() - t0

    print(f"  n={n_size}  ({2*n_size} flat modes)  "
          f"native(explicit)={t_native:.2f}s  "
          f"doubled(explicit)={t_double_explicit:.2f}s  "
          f"doubled(full,C++)={t_double_full:.2f}s")
