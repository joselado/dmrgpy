"""What a bare ARPACK call does to a degenerate level, and what dmrgpy does.

Ferromagnetic Heisenberg chain, 12 spin-1/2 sites -> dim 4096, which is
above `algebra.maxsize` (2000), so the ED path is iterative rather than
dense. With H = -sum S_i.S_{i+1} the ground level is exactly -2.75 and
exactly 13-fold degenerate (S=6), so asking for the lowest 8 levels must
give eight copies of the same energy.

A plain `eigsh(M, k=8)` does not: it stops once the Ritz pairs it holds have
converged, a condition one copy of a degenerate eigenvalue already
satisfies. That was dmrgpy's ED path, and it made `mode="ED"` look like it
returned distinct levels while `mode="DMRG"` returned multiplet members --
see docs/ed_vs_dmrg_degenerate_multiplets.md for how far that got.

`algebra.lowest_states` now deflates instead (see
`algebra._deflated_lowest_hermitian`), so the last two rows below match
dense. The bare-eigsh rows are kept precisely because they still fail: this
script's job is to show the difference, not to pass.
"""
import os
import sys

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
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
for label, w in [("dense eigvalsh", dense),
                 ("bare eigsh tol=0", t0),
                 ("bare eigsh tol=1e-7", t7),
                 ("algebra.lowest_states", ls),
                 ("get_excited(mode='ED')", ex)]:
    print("  %-26s copies of E0: %d/8   %s" % (
        label, copies(w), np.round(np.asarray(w).real, 6)))
