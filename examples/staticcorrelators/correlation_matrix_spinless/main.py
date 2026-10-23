# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Correlated orbitals of a spinless fermionic chain with nearest neighbor
# interaction U. We diagonalize the one-body correlation matrix
# <c^dagger_i c_j> of the ground state. Its eigenvectors are the natural
# orbitals, and eigenvalues away from 0 and 1 measure the correlations.
import numpy as np
from dmrgpy import fermionchain

n = 40 # number of sites
fc = fermionchain.Fermionic_Chain(n)
U = -0.95 # nearest neighbor interaction
h = 0
for i in range(n-1):
    h = h + fc.Cdag[i]*fc.C[i+1] # hopping
    h = h + U*(fc.N[i]-.5)*(fc.N[i+1]-.5) # interaction
h = h + h.get_dagger()
fc.maxm = 20 # bond dimension
fc.nsweeps = 10 # number of DMRG sweeps
fc.set_hamiltonian(h)
fc.get_gs()

out = fc.get_correlated_orbitals(ordered=False) # out[i][j], orbital i on site j

import matplotlib.pyplot as plt
plt.imshow(np.array(out), aspect="auto", origin="lower")
plt.colorbar(label="orbital weight")
plt.xlabel("site j")
plt.ylabel("orbital i")
plt.show()
