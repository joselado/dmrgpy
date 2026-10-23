# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Low energy excited states E_n of a Heisenberg chain with random
# couplings, H = sum_i J_i S_i . S_{i+1}. They give the excitation
# spectrum that inelastic neutron scattering or spin-polarized STM probes.
# We compute them with DMRG and compare with exact diagonalization (ED).
import numpy as np
from dmrgpy import spinchain

n = 8 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
np.random.seed(1) # same random couplings every time
h = 0
for i in range(n-1): # Heisenberg Hamiltonian with random couplings
    h = h + np.random.random()*sc.SS(i,i+1)
sc.set_hamiltonian(h)
sc.maxm = 40 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps

es1, ws1 = sc.get_excited_states(n=6, mode="DMRG") # 6 lowest states
es2, ws2 = sc.get_excited_states(n=6, mode="ED")
print("Max difference DMRG-ED =", np.max(np.abs(np.array(es1)-es2)))

import matplotlib.pyplot as plt
plt.plot(es2, "o", markersize=10, label="ED")
plt.plot(es1, "x", markersize=10, label="DMRG")
plt.xlabel("state index n")
plt.ylabel("energy E_n")
plt.legend()
plt.show()
