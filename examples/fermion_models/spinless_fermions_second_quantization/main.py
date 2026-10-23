# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Ground-state energy of a random spinless fermion chain written in second
# quantization, with first neighbor hopping, pairing and interaction,
# H = sum_i t_i c^+_i c_{i+1} + d_i c_i c_{i+1} + U_i n_i n_{i+1} + h.c.
# We compute it with DMRG and with exact diagonalization (ED).
import numpy as np
from dmrgpy import fermionchain

nfs = [4, 6, 8, 10] # numbers of sites
edmrg, eed = [], []
for n in nfs:
    np.random.seed(1) # same random numbers in every run
    fc = fermionchain.Fermionic_Chain(n)
    h = 0
    for i in range(n-1):
        h = h + fc.Cdag[i]*fc.C[i+1]*np.random.random() # hopping
        h = h + fc.C[i]*fc.C[i+1]*np.random.random() # pairing
        h = h + fc.N[i]*fc.N[i+1]*np.random.random() # interaction
    h = h + h.get_dagger() # make it Hermitian
    fc.set_hamiltonian(h)
    fc.maxm = 40 # bond dimension
    fc.nsweeps = 10 # number of DMRG sweeps
    edmrg.append(fc.gs_energy(mode="DMRG"))
    eed.append(fc.gs_energy(mode="ED"))
print("Max difference DMRG-ED =", np.max(np.abs(np.array(edmrg)-eed)))

import matplotlib.pyplot as plt
plt.plot(nfs, eed, linewidth=4, label="ED")
plt.plot(nfs, edmrg, "--", linewidth=2, label="DMRG")
plt.xlabel("number of sites")
plt.ylabel("ground-state energy")
plt.legend()
plt.show()
