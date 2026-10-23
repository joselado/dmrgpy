# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Ground state energy of the S=1/2 Heisenberg chain for chains too large
# for exact diagonalization. We compute it with the two ITensor backends,
# v2 and v3, and check that both give the same energy.
import numpy as np
from dmrgpy import spinchain

def energy(version, n):
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)], itensor_version=version)
    h = 0
    for i in range(n-1): # Heisenberg Hamiltonian
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = 40 # bond dimension
    sc.nsweeps = 10 # number of DMRG sweeps
    return sc.gs_energy(mode="DMRG")

ns = [16, 24, 32] # chain lengths
e2 = [energy(2, n) for n in ns]
e3 = [energy(3, n) for n in ns]
diff = np.max(np.abs(np.array(e2)-e3))
print("Max difference v2-v3 =", diff)
assert diff < 1e-3

import matplotlib.pyplot as plt
plt.plot(ns, e2, linewidth=4, label="v2")
plt.plot(ns, e3, "--", linewidth=2, label="v3")
plt.xlabel("chain length n")
plt.ylabel("ground state energy")
plt.legend()
plt.show()
