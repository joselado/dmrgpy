# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Fermionic parity of the ground state of a spinful chain with hopping,
# Hubbard U, s-wave pairing and a magnetic impurity on the first site.
# The parity is +1 (-1) for an even (odd) number of electrons. It is
# computed with DMRG, as a function of the length of the chain L.
import numpy as np
from dmrgpy import fermionchain

def parity(L, mu=0., U=0.5, Bz=1.0):
    fc = fermionchain.Spinful_Fermionic_Chain(L)
    h = 0
    for i in range(L-1): # hopping
        h = h + fc.Cdagup[i]*fc.Cup[i+1] + fc.Cdagdn[i]*fc.Cdn[i+1]
    for i in range(L): # chemical potential and pairing
        h = h + mu*fc.Cdagup[i]*fc.Cup[i] + mu*fc.Cdagdn[i]*fc.Cdn[i]
        h = h + mu*fc.Cup[i]*fc.Cdn[i]
    for i in range(L-1): # Hubbard interaction
        h = h + U*(fc.Nup[i]-.5)*(fc.Ndn[i]-.5)
    h = h + Bz*fc.Sz[0] # magnetic impurity
    fc.set_hamiltonian(h + h.get_dagger())
    fc.maxm = 40 # bond dimension
    fc.nsweeps = 20 # number of DMRG sweeps
    wf = fc.get_gs(mode="DMRG")
    return wf.get_fermionic_parity(fpmode="iterative").real

Ls = [4, 8, 12, 16, 20] # lengths of the chain
ps = [parity(L) for L in Ls]
print("Parities =", ps)

import matplotlib.pyplot as plt
plt.plot(Ls, ps, "o")
plt.xlabel("chain length L")
plt.ylabel("fermionic parity")
plt.show()
