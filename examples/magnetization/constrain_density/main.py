# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Fixed density with a Lagrange multiplier. We add a penalty
# lambda*(N - N0)^2 to a random hopping Hamiltonian, so that the ground
# state has a total number of electrons N close to the target N0. This is
# a simple way to work at fixed density without a conserved sector.
import numpy as np
from dmrgpy import fermionchain

n = 5 # number of sites
fc = fermionchain.Spinful_Fermionic_Chain(n)
m = np.random.random((n, n)) + 1j*np.random.random((n, n))
m = m + np.conjugate(m).T # random Hermitian hopping matrix
h0 = 0
for i in range(n):
    for j in range(n):
        h0 = h0 + m[i,j]*fc.Cdag[i]*fc.C[j]

lamb = 100 # Lagrange multiplier for the density
targets = range(0, 2*n+1) # target number of electrons
densities = []
for N0 in targets:
    Nop = sum(fc.N) - N0
    fc.set_hamiltonian(h0 + lamb*Nop*Nop)
    wf = fc.get_gs(mode="ED") # ground state with exact diagonalization
    densities.append(wf.dot(sum(fc.N)*wf).real)
print("Max difference =", np.max(np.abs(np.array(densities)-list(targets))))

import matplotlib.pyplot as plt
plt.plot(targets, targets, linewidth=4, label="target")
plt.plot(targets, densities, "--", linewidth=2, label="ground state")
plt.xlabel("target number of electrons N0")
plt.ylabel("<N>")
plt.legend()
plt.show()
