# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Leaving a conserved sector without losing the state. We solve the spinless
# t-V chain at half filling with the particle number fixed (fast and exact in
# the sector), and then call promote_to_dense(). Now an operator that changes
# the particle number, like c_i, can be applied. The weight |c_i|GS>|^2 is
# the occupation <n_i>, the charge density wave that photoemission can see.
import numpy as np
from dmrgpy import fermionchain

n = 10 # number of sites
V = 2.0 # nearest neighbor repulsion
fc = fermionchain.Fermionic_Chain(n)
h = 0
for i in range(n-1):
    h = h + fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i] + V*fc.N[i]*fc.N[i+1]
fc.set_hamiltonian(h)
fc.maxm = 60 # bond dimension
fc.nsweeps = 10 # number of DMRG sweeps

fc.set_conserved_sector(Nf=n//2) # fix the number of particles
print("Energy in the sector =", fc.gs_energy())
fc.promote_to_dense() # leave the sector, keep the state
wf = fc.wf0.copy()

# c_i changes the particle number, it is allowed only after promoting
weights = np.array([fc.applyoperator(fc.C[i], wf).norm()**2 for i in range(n)])
occ = np.array([fc.vev(fc.N[i]).real for i in range(n)])
print("Max difference |c_i GS|^2 - <n_i> =", np.max(np.abs(weights-occ)))

import matplotlib.pyplot as plt
plt.plot(range(n), occ, "o", markersize=12, label="<n_i>")
plt.plot(range(n), weights, "x", markersize=10, label="|c_i GS|^2")
plt.xlabel("site i")
plt.ylabel("occupation")
plt.legend()
plt.show()
