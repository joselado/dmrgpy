# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Hopping correlator <Cdag_0 C_j> in the ground state of a chain of
# spinless fermions, H = sum_i Cdag_i C_{i+1} + h.c. Fermionic operators
# on different sites carry a Jordan-Wigner string, and this example checks
# that DMRG gives the same correlator as exact diagonalization (ED).
import numpy as np
from dmrgpy import fermionchain

n = 8 # number of sites
fc = fermionchain.Fermionic_Chain(n)

# tight-binding Hamiltonian with first-neighbor hopping
h = 0
for i in range(n-1):
    h = h + fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i]
fc.set_hamiltonian(h)
fc.maxm = 40 # bond dimension
fc.nsweeps = 10 # number of DMRG sweeps

# correlation between the first site and site j
cdmrg = np.array([fc.vev(fc.Cdag[0]*fc.C[j], mode="DMRG") for j in range(n)])
ced = np.array([fc.vev(fc.Cdag[0]*fc.C[j], mode="ED") for j in range(n)])
err = np.max(np.abs(cdmrg - ced))
print("Max difference DMRG-ED =", err)
assert err < 1e-3

import matplotlib.pyplot as plt
plt.plot(range(n), ced.real, "o", markersize=10, label="ED")
plt.plot(range(n), cdmrg.real, "x", markersize=10, label="DMRG")
plt.xlabel("site j")
plt.ylabel("<Cdag_0 C_j>")
plt.legend()
plt.show()
