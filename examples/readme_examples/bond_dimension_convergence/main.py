# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Convergence of DMRG with the bond dimension, for the S=1/2 Heisenberg
# chain. We plot the energy per site, and the energy fluctuation
# sqrt(<H^2>-<H>^2) per site as an error bar. The fluctuation is zero
# for an exact eigenstate, so it tells how good the state is.
import numpy as np
from dmrgpy import spinchain

n = 30 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
h = 0
for i in range(n-1): # Heisenberg Hamiltonian
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]

bds = range(3,20,2) # bond dimensions
es, des = [], []
for maxm in bds:
    sc.set_hamiltonian(h)
    sc.maxm = maxm # bond dimension
    sc.nsweeps = 10 # number of DMRG sweeps
    e = sc.gs_energy()
    wf = sc.get_gs()
    de = np.sqrt(np.abs(wf.dot(h*(h*wf)) - e**2)) # energy fluctuation
    es.append(e/n)
    des.append(de/n)
    print(maxm, e/n, de/n)

import matplotlib.pyplot as plt
plt.errorbar(bds, es, yerr=des, marker="o", capsize=5)
plt.xlabel("bond dimension")
plt.ylabel("energy per site")
plt.show()
