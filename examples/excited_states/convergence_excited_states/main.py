# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Excited states of the S=1/2 Heisenberg chain with DMRG. Each excited
# state is found with a penalty on the lower states. We check how the
# error of the first six energies, compared with ED, decreases with the
# number of DMRG sweeps.
import numpy as np
from dmrgpy import spinchain

n = 6 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 40 # bond dimension
e_ed = sc.get_excited(n=6, mode="ED") # exact energies

nsweeps = [4, 6, 8, 10, 20, 30]
errors = []
for ns in nsweeps:
    sc.nsweeps = ns # number of DMRG sweeps
    e = sc.get_excited(n=6, scale=10) # scale is the penalty of the lower states
    errors.append(np.mean(np.abs(np.array(e_ed) - np.array(e))))
    print(ns, errors[-1])

import matplotlib.pyplot as plt
plt.plot(nsweeps, errors, "o-")
plt.yscale("log")
plt.xlabel("Number of sweeps")
plt.ylabel("Mean error of the energies")
plt.show()
