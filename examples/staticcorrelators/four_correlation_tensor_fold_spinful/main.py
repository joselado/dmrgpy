"""ctmode="fold" for the 4-point correlator tensor on native spinful sites.

Spinful_Fermionic_Chain_Native (ITensor ElectronSite) has no "sweep": its
Cup/Cdn/Cdagup/Cdagdn operators are not the single unified fermion flavour
that path assumes. Until "fold" existed the only options were a per-tuple
AutoMPO build -- ctmode="explicit" on the pure-Python backend, ctmode="full"
on the C++ one -- which spends essentially all its time compiling an MPO over
the whole chain and sweeping it, for an operator that touches at most 4 sites.

"fold" evaluates each tuple as a local operator fold instead. This script
checks it against those exact per-tuple oracles and plots how the saving
grows with system size, on whichever backends are available.

Run: python3 main.py    (add MKL_NUM_THREADS=1 for stable timings -- see
docs/user_guide.md's BLAS-threads section)
"""
import os
import sys
import time

sys.path.append(os.path.join(os.path.dirname(__file__), "..", "..", "..", "src"))

import numpy as np

from dmrgpy import cppext, fermionchain


def build(n, itensor_version):
    """Hubbard-like spinful chain: nearest-neighbour hopping per spin, on-site
    U, plus a small potential on mode 0 to break particle-hole symmetry (so
    the tensor has non-trivial structure rather than mostly zeros)."""
    fc = fermionchain.Spinful_Fermionic_Chain_Native(
        n, itensor_version=itensor_version)
    h = 0
    for i in range(n - 1):
        for f in (0, 1):                       # 0 = up, 1 = down
            a, b = 2 * i + f, 2 * (i + 1) + f
            h = h + fc.Cdag[a] * fc.C[b] + fc.Cdag[b] * fc.C[a]
    for i in range(n):                         # U n_up n_dn
        h = h + 2.0 * (fc.Cdag[2*i] * fc.C[2*i]) * (fc.Cdag[2*i+1] * fc.C[2*i+1])
    h = h + 0.3 * (fc.Cdag[0] * fc.C[0])
    fc.set_hamiltonian(h)
    fc.maxm = 20
    fc.gs_energy()
    return fc


backends = [("python", "explicit")]
if cppext.available(3):
    backends.append((3, "full"))

results = {}
for version, slow_mode in backends:
    label = f'itensor_version={version!r}'
    print(f"\n### {label}: fold vs {slow_mode} ###")
    ns, t_fold, t_slow = [], [], []
    for n in (2, 3, 4, 5):
        if version == 3 and n < 3:
            continue          # ITensor v3's two-site DMRG needs >= 3 sites
        wf = build(n, version).get_gs()

        t0 = time.time()
        fold = wf.get_four_correlation_tensor(ctmode="fold")
        dt_fold = time.time() - t0

        t0 = time.time()
        slow = wf.get_four_correlation_tensor(ctmode=slow_mode)
        dt_slow = time.time() - t0

        err = np.abs(fold - slow).max()
        print(f"  sites={n} modes={2*n:2d}: max|fold-{slow_mode}| = {err:.2e}   "
              f"fold {dt_fold:6.2f}s   {slow_mode} {dt_slow:6.2f}s   "
              f"({dt_slow/max(dt_fold, 1e-9):4.1f}x)")
        assert err < 1e-10, (n, err)   # both are the same tensor, exactly
        ns.append(2 * n)
        t_fold.append(dt_fold)
        t_slow.append(dt_slow)
    results[label] = (ns, t_fold, t_slow, slow_mode)

import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, len(results), figsize=(5.5 * len(results), 4),
                          squeeze=False)
for ax, (label, (ns, t_fold, t_slow, slow_mode)) in zip(axes[0], results.items()):
    ax.plot(ns, t_slow, marker="o", label=f'ctmode="{slow_mode}"')
    ax.plot(ns, t_fold, marker="s", label='ctmode="fold"')
    ax.set_yscale("log")
    ax.set_xlabel("fermionic modes (2 x sites)")
    ax.set_ylabel("time [s]")
    ax.set_title(label)
    ax.legend()
plt.tight_layout()
plt.show()
