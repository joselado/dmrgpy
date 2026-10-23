# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Lowest energies of a non-Hermitian fermionic chain, a tight-binding
# chain with a staggered imaginary onsite term, H = sum_i c^dag_i c_{i+1}
# + h.c. + sum_i (-1)^i 0.6j n_i. The energies are complex, and they
# describe the lifetime (imaginary part) of the excitations. We compute
# them with DMRG and compare with exact diagonalization (ED).
import numpy as np
from dmrgpy import fermionchain

n = 6 # number of sites
fc = fermionchain.Fermionic_Chain(n)
h = 0
for i in range(n-1):
    h = h + fc.Cdag[i]*fc.C[i+1]
h = h + h.get_dagger() # Hermitian hopping
for i in range(n):
    h = h + (-1)**i*0.6j*fc.Cdag[i]*fc.C[i] # staggered imaginary onsite term
fc.set_hamiltonian(h)
fc.maxm = 40 # bond dimension
fc.nsweeps = 10 # number of DMRG sweeps

e1 = np.array(fc.get_excited(mode="DMRG", n=4))
e2 = np.array(fc.get_excited(mode="ED", n=4))
print("Energies ED =", np.round(e2, 3))

import matplotlib.pyplot as plt
plt.scatter(e2.real, e2.imag, s=150, marker="o", label="ED")
plt.scatter(e1.real, e1.imag, s=60, marker="x", label="DMRG")
plt.xlabel("Re(E)")
plt.ylabel("Im(E)")
plt.legend()
plt.show()
