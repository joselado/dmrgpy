# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Energy fluctuation per site of the ground state of the S=1 Heisenberg
# chain, ||(H-<H>)|psi>||/n, which is zero for an exact eigenstate.
# With gs_energy(maxde=...) the bond dimension is increased until
# the fluctuation per site is below the value requested.
import numpy as np
from dmrgpy import spinchain

n = 10 # number of sites
sc = spinchain.Spin_Chain(["S=1" for i in range(n)])

# Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.maxm = 10 # starting bond dimension
sc.nsweeps = 15 # number of DMRG sweeps

maxdes = [1e-1, 1e-3, 1e-6] # requested fluctuation per site
fluctuations = []
for maxde in maxdes:
    sc.set_hamiltonian(h) # start again for each maxde
    sc.gs_energy(maxde=maxde)
    fluctuations.append(sc.gs_energy_fluctuation()/n)
    print("maxde", maxde, "reached", fluctuations[-1])

import matplotlib.pyplot as plt
plt.loglog(maxdes, fluctuations, "o-", label="reached")
plt.loglog(maxdes, maxdes, "k--", label="requested")
plt.xlabel("requested fluctuation per site")
plt.ylabel("fluctuation per site")
plt.legend()
plt.show()
