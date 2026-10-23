# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Expectation values <O> of the parafermion operators Chi and Psi in the
# ground state of a random Z3 parafermion chain. Parafermions generalize
# Majorana fermions, and their operators fulfill a Z3 algebra. We compute
# the expectation values with DMRG and compare with exact diagonalization.
import numpy as np
from dmrgpy import parafermionchain

n = 3 # number of sites
np.random.seed(1) # same random Hamiltonian in every run
pc = parafermionchain.Parafermionic_Chain(n)
h = 0
for i in range(n):
    for j in range(n):
        h = h + pc.Sig[i]*pc.Tau[j]*np.random.random()
        h = h + 1j*pc.Sig[i]*pc.Sigd[j]*np.random.random()
        h = h + pc.Tau[i]*pc.Tau[j]*np.random.random()
        h = h + 1j*pc.Tau[i]*pc.Taud[j]*np.random.random()
h = h + h.get_dagger() # make it Hermitian
pc.set_hamiltonian(h)
pc.maxm = 60 # bond dimension
pc.nsweeps = 40 # number of DMRG sweeps

ops = pc.Chi + pc.Psi # all the parafermion operators
c1 = np.array([pc.vev(o, mode="ED") for o in ops])
c2 = np.array([pc.vev(o, mode="DMRG") for o in ops])
print("Max difference DMRG-ED =", np.max(np.abs(c1-c2)))

import matplotlib.pyplot as plt
plt.plot(range(len(ops)), c1.real, "o", markersize=10, label="ED")
plt.plot(range(len(ops)), c2.real, "x", markersize=10, label="DMRG")
plt.xlabel("operator")
plt.ylabel("<O>")
plt.legend()
plt.show()
