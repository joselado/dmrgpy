# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Ground state energy of the anisotropic Heisenberg chain,
# H = sum_i Sz_i Sz_{i+1} + alpha (Sx_i Sx_{i+1} + Sy_i Sy_{i+1}).
# alpha=0 is the Ising chain and alpha=1 is the Heisenberg chain.
# We compute the energy with DMRG and check it against exact diagonalization (ED).
import numpy as np
from dmrgpy import spinchain

n = 6 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
sc.maxm = 20 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps

alphas = np.linspace(0, 1, 6)
e_dmrg, e_ed = [], []
for alpha in alphas:
    h = 0
    for i in range(n-1):
        h = h + sc.Sz[i]*sc.Sz[i+1]
        h = h + alpha*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1])
    sc.set_hamiltonian(h)
    e_dmrg.append(sc.gs_energy(mode="DMRG"))
    e_ed.append(sc.gs_energy(mode="ED"))
print("Max difference DMRG-ED =", np.max(np.abs(np.array(e_dmrg)-e_ed)))

import matplotlib.pyplot as plt
plt.plot(alphas, e_ed, linewidth=4, label="ED")
plt.plot(alphas, e_dmrg, "--", linewidth=2, label="DMRG")
plt.xlabel("alpha")
plt.ylabel("Ground state energy")
plt.legend()
plt.show()
