# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Spinful fermion chain with hopping, Hubbard interaction U and
# nearest-neighbor repulsion V (the UV chain), all in units of the hopping t=1.
# We compute the ground state energy for several lengths with DMRG
# and compare with exact diagonalization (ED).
import numpy as np
from dmrgpy import fermionchain

def energies(n): # DMRG and ED ground state energies
    fc = fermionchain.Spinful_Fermionic_Chain(n)
    h = 0
    for i in range(n-1): # hopping
        h = h + fc.Cdagup[i]*fc.Cup[i+1] + fc.Cdagdn[i]*fc.Cdn[i+1]
    for i in range(n): # Hubbard interaction, U=1
        h = h + (fc.Nup[i]-.5)*(fc.Ndn[i]-.5)
    for i in range(n-1): # nearest-neighbor repulsion, V=1
        h = h + (fc.Nup[i]+fc.Ndn[i])*(fc.Nup[i+1]+fc.Ndn[i+1])
    h = h + h.get_dagger()
    fc.maxm = 40 # bond dimension
    fc.nsweeps = 40 # number of DMRG sweeps
    fc.set_hamiltonian(h)
    return fc.gs_energy(mode="DMRG"), fc.gs_energy(mode="ED")

ns = [2, 3, 4, 5] # number of spinful sites
e = np.array([energies(n) for n in ns])
print("Max difference DMRG-ED =", np.max(np.abs(e[:,0]-e[:,1])))

import matplotlib.pyplot as plt
plt.plot(ns, e[:,1], linewidth=4, label="ED")
plt.plot(ns, e[:,0], "--", linewidth=2, label="DMRG")
plt.xlabel("number of sites")
plt.ylabel("ground state energy")
plt.legend()
plt.show()
