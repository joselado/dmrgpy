# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Density <Cdag_i C_i> of a spinless fermionic chain with random complex
# hopping between all sites, built from MultiOperator objects. We compute
# it with DMRG and with exact diagonalization (ED), and compare both.
import numpy as np
from dmrgpy import fermionchain
from dmrgpy import multioperator

n = 6 # number of sites
np.random.seed(1) # same random hopping every run
fc = fermionchain.Fermionic_Chain(n)
fc.maxm = 40 # bond dimension
fc.nsweeps = 60 # number of DMRG sweeps
m = np.random.random((n,n)) + 1j*np.random.random((n,n))
m = m + m.conj().T # Hermitian hopping matrix
h = 0
for i in range(n):
    for j in range(n):
        h = h + m[i,j]*fc.Cdag[i]*fc.C[j]
fc.set_hamiltonian(h)

den = [multioperator.obj2MO([["Cdag",i]])*multioperator.obj2MO([["C",i]]) for i in range(n)]
d1 = [fc.vev(d, mode="DMRG").real for d in den]
d2 = [fc.vev(d, mode="ED").real for d in den]
print("Max difference DMRG-ED =", np.max(np.abs(np.array(d1)-d2)))

import matplotlib.pyplot as plt
plt.plot(range(n), d2, linewidth=4, label="ED")
plt.plot(range(n), d1, "--", linewidth=2, label="DMRG")
plt.xlabel("site i")
plt.ylabel("<Cdag_i C_i>")
plt.legend()
plt.show()
