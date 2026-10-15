# Energy per site of the infinite spin-1/2 Heisenberg chain against the bond
# dimension (iDMRG), and the spin correlator at long distance.
import numpy as np
from dmrgpy import infinitechain

maxms = [8, 12, 16, 24, 32, 48, 64, 96]
rmax = 40
e_exact = 0.25 - np.log(2)  # Bethe ansatz


def heisenberg(cell, m):
    ic = infinitechain.Infinite_Spin_Chain(cell)
    h = 0
    for o in ("Sx", "Sy", "Sz"):
        oc, orr = getattr(ic, o + "C"), getattr(ic, o + "R")
        if len(cell) == 1:
            h = h + oc[0]*orr[0]
        else:
            h = h + oc[0]*oc[1] + oc[1]*orr[0]
    ic.set_hamiltonian(h)
    ic.gs_method, ic.maxm, ic.maxiter, ic.etol = "idmrg", m, 400, 1e-10
    return ic


energy, converged = [], []
for m in maxms:
    ic = heisenberg(["1/2"], m)
    energy.append(ic.gs_energy())
    converged.append(bool(ic.converged))

# At finite maxm the iDMRG state of this critical chain alternates between
# strong and weak bonds, so the correlator is taken on a two-site cell and
# averaged over its two sites.
ic2 = heisenberg(["1/2", "1/2"], maxms[-1])
e_uc2 = ic2.gs_energy()
r = np.arange(1, rmax + 1)
czz = np.mean([[ic2.correlator("Sz", i, "Sz", int(x)).real for x in r]
               for i in (0, 1)], axis=0)

np.savez("data.npz", maxm=maxms, energy=energy, converged=converged,
         e_exact=e_exact, e_uc2=e_uc2, r=r, czz=czz)
