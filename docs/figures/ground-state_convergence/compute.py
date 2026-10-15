# Energy error and energy fluctuation of the spin-1/2 Heisenberg chain
# against the bond dimension, with exact diagonalization as the reference.
import numpy as np
from dmrgpy import spinchain

n = 16
maxms = [2, 4, 8, 12, 16, 24, 32, 48, 64, 96]


def heisenberg():
    sc = spinchain.Spin_Chain(["S=1/2"]*n)
    h = sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
            for i in range(n-1))
    sc.set_hamiltonian(h)
    return sc


e_ed = heisenberg().gs_energy(mode="ED")
energy, de = [], []
for maxm in maxms:
    sc = heisenberg()  # a fresh chain for each bond dimension
    sc.maxm, sc.nsweeps, sc.cutoff = maxm, 20, 1e-14
    energy.append(sc.gs_energy())
    de.append(sc.gs_energy_fluctuation())

np.savez("data.npz", n=n, maxm=maxms, energy=energy, de=de, e_ed=e_ed)
