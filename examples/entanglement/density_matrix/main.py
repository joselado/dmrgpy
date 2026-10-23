# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Eigenvalues of the reduced density matrix of the full chain, as a function
# of the interaction U, for spinless fermions at half filling,
# H = sum_i c^+_i c_{i+1} + h.c. + U (n_i - 1/2)(n_{i+1} - 1/2).
# For a pure state of the whole chain, the eigenvalues show how mixed
# the state is, and so how entangled it is with the part that is traced out.
import numpy as np
import scipy.linalg as lg
from dmrgpy import fermionchain

n = 6 # number of sites
fc = fermionchain.Fermionic_Chain(n)

def spectrum(U):
    h = 0
    for i in range(n-1):
        h = h + fc.Cdag[i]*fc.C[i+1]
        h = h + U*(fc.N[i]-0.5)*(fc.N[i+1]-0.5)
    h = h + h.get_dagger()
    fc.set_hamiltonian(h)
    wf = fc.get_gs(mode="ED") # ground state
    dm = wf.get_dm(inds=range(n)) # density matrix
    return lg.eigvalsh(dm)

import matplotlib.pyplot as plt
for U in np.linspace(0.0, 20.0, 20):
    y = spectrum(U)
    plt.scatter(0*y+U, y, c="black")
plt.xlabel("interaction U")
plt.ylabel("eigenvalues of the density matrix")
plt.show()
