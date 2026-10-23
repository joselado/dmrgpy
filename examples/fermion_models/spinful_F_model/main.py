# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Chain with two kinds of fermions on each site: spinful electrons with
# a Hubbard interaction, and one spinless pseudofermion F. The F fermion
# hops along the chain and hybridizes with the electron density. We
# compute the ground state energy as a function of the chain length.
import numpy as np
from dmrgpy import fermionchain

def get_energy(n):
    fc = fermionchain.Spinful_F_Fermionic_Chain(n)
    h = 0
    for i in range(n-1): # hopping of the electrons and of F
        h = h + fc.Cdagup[i]*fc.Cup[i+1] + fc.Cdagdn[i]*fc.Cdn[i+1]
        h = h + fc.Fdag[i]*fc.F[i+1]
    for i in range(n):
        h = h + (fc.Nup[i]-.5)*(fc.Ndn[i]-.5) # Hubbard interaction
        h = h + fc.Fdag[i]*(fc.Nup[i]+fc.Ndn[i]) # hybridization
    fc.set_hamiltonian(h + h.get_dagger())
    return fc.gs_energy()

ns = [2, 3, 4, 5] # number of sites
es = [get_energy(n) for n in ns]
print("Energies =", es)

import matplotlib.pyplot as plt
plt.plot(ns, es, "o-")
plt.xlabel("number of sites")
plt.ylabel("ground state energy")
plt.show()
