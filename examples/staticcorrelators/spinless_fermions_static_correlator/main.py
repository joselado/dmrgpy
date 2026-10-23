# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Static correlators <c^+_i c_j> of spinless fermions with a random
# hopping matrix. They give the density (i=j) and the bond order (i!=j)
# of the ground state. We compute them with DMRG and with exact
# diagonalization (ED), and compare both.
import numpy as np
from dmrgpy import fermionchain

n = 4 # number of orbitals
fc = fermionchain.Fermionic_Chain(n)
np.random.seed(3)
h = 0
for i in range(n):
    for j in range(n): # random hopping
        h = h + np.random.random()*fc.Cdag[i]*fc.C[j]
h = h + h.get_dagger()
fc.set_hamiltonian(h)
fc.maxm = 40 # bond dimension
fc.nsweeps = 10 # number of DMRG sweeps

pairs = [(i, j) for i in range(n) for j in range(i, n)]
c1 = np.array([fc.vev(fc.Cdag[i]*fc.C[j], mode="DMRG") for i, j in pairs])
c2 = np.array([fc.vev(fc.Cdag[i]*fc.C[j], mode="ED") for i, j in pairs])
print("Max difference DMRG-ED =", np.max(np.abs(c1-c2)))

import matplotlib.pyplot as plt
plt.plot(range(len(pairs)), c2.real, "o", markersize=10, label="ED")
plt.plot(range(len(pairs)), c1.real, "x", markersize=10, label="DMRG")
plt.xlabel("pair number (i,j)")
plt.ylabel("<c^+_i c_j>")
plt.legend()
plt.show()
