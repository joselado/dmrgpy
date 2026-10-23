# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Ground state energy of the S=1/2 Heisenberg chain with the two C++
# DMRG backends, ITensor v2 and v3, as a function of the chain length.
# Both should give the same energy for every size.
import numpy as np
from dmrgpy import spinchain

def get_energy(version, n):
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)], itensor_version=version)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = 40 # bond dimension
    sc.nsweeps = 10 # number of DMRG sweeps
    return sc.gs_energy()

ns = [4, 6, 8, 10, 12] # chain lengths
e2 = [get_energy(2, n) for n in ns]
e3 = [get_energy(3, n) for n in ns]
print("Max difference v2-v3 =", np.max(np.abs(np.array(e2)-e3)))

import matplotlib.pyplot as plt
plt.plot(ns, e2, linewidth=4, label="ITensor v2")
plt.plot(ns, e3, "--", linewidth=2, label="ITensor v3")
plt.xlabel("number of sites")
plt.ylabel("ground state energy")
plt.legend()
plt.show()
