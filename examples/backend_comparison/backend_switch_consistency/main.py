# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Switching the backend of a chain. We build a Heisenberg chain, compute
# the ground state energy with ITensor v2, switch the same chain to v3
# without setting the Hamiltonian again, and compare with exact
# diagonalization (ED). The energy after the switch must be a new v3
# result and not the old v2 one.
import numpy as np
from dmrgpy import spinchain

ns = [3,4,5,6] # chain lengths
e2, e3, eed = [], [], []
for n in ns:
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)], itensor_version=2)
    h = 0
    for i in range(n-1): # Heisenberg Hamiltonian
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = 20 # bond dimension
    sc.nsweeps = 10 # number of DMRG sweeps
    e2.append(sc.gs_energy(mode="DMRG")) # energy with v2
    sc.setup_cpp(3) # switch the same chain to v3
    e3.append(sc.gs_energy(mode="DMRG")) # energy with v3
    eed.append(sc.gs_energy(mode="ED")) # exact energy
diff = np.max(np.abs(np.array(e3)-eed))
print("Max difference v3-ED =", diff)
assert diff < 1e-3

import matplotlib.pyplot as plt
plt.plot(ns, eed, linewidth=4, label="ED")
plt.plot(ns, e2, "o", label="v2")
plt.plot(ns, e3, "x--", label="v3 after the switch")
plt.xlabel("chain length n")
plt.ylabel("ground state energy")
plt.legend()
plt.show()
