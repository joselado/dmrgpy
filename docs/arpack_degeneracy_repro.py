"""Minimal reproducer: dmrgpy's mode="ED" under-reports degenerate levels
above algebra.maxsize, because algebra.lowest_states falls back to ARPACK
there and Lanczos from a single start vector cannot resolve multiplicity.

Ferromagnetic Heisenberg chain, 12 spin-1/2 sites -> dim 4096 (above
maxsize=2000) with a 13-fold degenerate ground multiplet (S=6). Asking for
the lowest 8 levels should give eight copies of the same energy.

Note the tolerance dependence: dmrgpy passes tol=1e-7 to eigsh, which makes
the collapse markedly worse than ARPACK's own default (tol=0).
"""
import sys
sys.path.insert(0, "/home/joselado/Documents/programs/dmrgpy-idmrg/src")
import numpy as np
import scipy.sparse.linalg as sl
from dmrgpy import spinchain
from dmrgpy.algebra import algebra

n = 12
sc = spinchain.Spin_Chain(["1/2"] * n)
h = 0
for i in range(n - 1):
    h = h - (sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
sc.set_hamiltonian(h)
M = sc.get_ED_obj().get_hamiltonian()

def copies(w, e0=-2.75):
    return int(np.sum(np.abs(np.asarray(w).real - e0) < 1e-8))

dense = np.sort(np.linalg.eigvalsh(np.array(M.todense())))[:8]
t0 = np.sort(sl.eigsh(M, k=8, which="SA", return_eigenvectors=False))
t7 = np.sort(sl.eigsh(M, k=8, which="SA", tol=1e-7, maxiter=int(1e6),
                       return_eigenvectors=False))
ls, _ = algebra.lowest_states(M, n=8)
ex = np.array(sc.get_excited(n=8, mode="ED"))

print("dim %d, ground multiplet is 13-fold; asking for the lowest 8:" % M.shape[0])
for label, w in [("dense eigvalsh", dense), ("eigsh tol=0", t0),
                 ("eigsh tol=1e-7 (dmrgpy's)", t7),
                 ("algebra.lowest_states", ls), ("get_excited(mode='ED')", ex)]:
    print("  %-26s copies of E0: %d/8   %s" % (label, copies(w), np.round(np.asarray(w).real, 6)))
