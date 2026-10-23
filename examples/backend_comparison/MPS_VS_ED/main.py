# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Ground state energy of random spin chains, with random spin sizes
# (S=1 to S=3/2 here) and random exchange between all pairs of sites.
# The energy from DMRG (matrix product states) is compared with the
# exact one from exact diagonalization (ED).
import numpy as np
from dmrgpy import spinchain

np.random.seed(1)
n = 5 # number of sites
e_dmrg, e_ed = [], []
for k in range(4): # four random chains
    spins = [np.random.randint(2, 5) for i in range(n)] # 2S+1 of each site
    sc = spinchain.Spin_Chain(spins)
    h = 0
    for i in range(n):
        for j in range(n): # random exchange between all pairs
            h = h + np.random.random()*(sc.Sx[i]*sc.Sx[j] + sc.Sy[i]*sc.Sy[j] + sc.Sz[i]*sc.Sz[j])
    sc.set_hamiltonian(h)
    sc.maxm = 60 # bond dimension
    sc.nsweeps = 10 # number of DMRG sweeps
    e_dmrg.append(sc.gs_energy(mode="DMRG"))
    e_ed.append(sc.gs_energy(mode="ED"))
print("Max difference DMRG-ED =", np.max(np.abs(np.array(e_dmrg)-e_ed)))

import matplotlib.pyplot as plt
plt.plot(e_ed, e_ed, linewidth=4, label="ED")
plt.plot(e_ed, e_dmrg, "x", markersize=10, label="DMRG")
plt.xlabel("energy ED")
plt.ylabel("energy")
plt.legend()
plt.show()
