# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Lowest excited states of the Heisenberg chain with DMRG. Each state is
# found by a sweep with the previous ones projected out, so errors pile
# up. The purification step cleans them. We show the mean error of the
# energies with respect to ED, with and without it, against the sweeps.
import numpy as np
from dmrgpy import spinchain

n = 8 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
h = 0
for i in range(n-1):
    h = h + sc.SS(i,i+1) # S_i . S_{i+1}
sc.set_hamiltonian(h)
sc.maxm = 10 # bond dimension (small on purpose)

eed = sc.get_excited(n=10, mode="ED") # reference energies
sweeps = [2, 4, 6] # numbers of DMRG sweeps
err, errpure = [], []
for nsweeps in sweeps:
    sc.nsweeps = nsweeps
    e = sc.get_excited(n=10, mode="DMRG", purify=False)
    err.append(np.mean(np.abs(e - eed)))
    e = sc.get_excited(n=10, mode="DMRG", purify=True)
    errpure.append(np.mean(np.abs(e - eed)))

import matplotlib.pyplot as plt
plt.semilogy(sweeps, err, "o-", label="no purification")
plt.semilogy(sweeps, errpure, "o-", label="purification")
plt.xlabel("number of sweeps")
plt.ylabel("mean error of the energies")
plt.legend()
plt.show()
