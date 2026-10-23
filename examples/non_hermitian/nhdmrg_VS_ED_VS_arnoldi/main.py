# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Non-Hermitian DMRG for an interacting fermionic chain with a staggered
# imaginary potential. The spectrum is complex, and the ground state is
# the eigenvalue with the smallest real part. We compute it with
# nhdmrg() and compare with the exact spectrum from ED.
import numpy as np
from dmrgpy import fermionchain

n = 6 # number of sites
fc = fermionchain.Fermionic_Chain(n)
h = 0
for i in range(n-1):
    h = h + fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i] # hopping
    h = h + (fc.N[i]-0.5)*(fc.N[i+1]-0.5) # interaction
for i in range(n):
    h = h + 1j*(-1)**i*fc.Cdag[i]*fc.C[i] # staggered imaginary potential
fc.set_hamiltonian(h)
fc.maxm = 40 # bond dimension
fc.nsweeps = 10 # number of DMRG sweeps

es = np.array(fc.get_excited(mode="ED", n=6)) # lowest exact energies
e, psil, psir = fc.nhdmrg() # energy, left and right eigenvectors
print("NH-DMRG energy:", e)
print("Distance to the closest exact energy =", np.min(np.abs(es-e)))
assert abs(e.real-es[0].real) < 1e-6

import matplotlib.pyplot as plt
plt.scatter(es.real, es.imag, s=100, label="ED")
plt.scatter([e.real], [e.imag], marker="x", s=100, label="NH-DMRG")
plt.xlabel("Re(E)")
plt.ylabel("Im(E)")
plt.legend()
plt.show()
