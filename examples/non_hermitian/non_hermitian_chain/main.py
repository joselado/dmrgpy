# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Low energy states of a non-Hermitian interacting fermionic chain,
# H = sum_ij M_ij c^dagger_i c_j + V sum_i (n_i-1/2)(n_{i+1}-1/2),
# with an imaginary staggered onsite term. The energies are complex. We
# find them with a non-Hermitian Arnoldi method on MPS and compare with ED.
import numpy as np
from dmrgpy import fermionchain, mpsalgebra

n = 4 # number of sites
fc = fermionchain.Fermionic_Chain(n)
h = 0
for i in range(n-1): # hopping and interaction
    h = h + fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i]
    h = h + (fc.N[i]-0.5)*(fc.N[i+1]-0.5)
for i in range(n): # imaginary staggered onsite term
    h = h + 1j*(-1)**i*fc.Cdag[i]*fc.C[i]
fc.set_hamiltonian(h)
fc.maxm = 40 # bond dimension

e_ed = np.array(fc.get_excited(mode="ED", n=3))
e_mps, wf = mpsalgebra.lowest_energy_non_hermitian_arnoldi(fc, h, n=3)
e_mps = np.array(e_mps)
print("ED energies", e_ed)
print("MPS energies", e_mps)

import matplotlib.pyplot as plt
plt.scatter(e_ed.real, e_ed.imag, s=200, label="ED")
plt.scatter(e_mps.real, e_mps.imag, s=60, label="MPS")
plt.xlabel("Re(E)")
plt.ylabel("Im(E)")
plt.legend()
plt.show()
