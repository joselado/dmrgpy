# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Ground state energy of the S=3/2 Heisenberg chain with three DMRG
# backends: ITensor v2, ITensor v3 and the pure Python one. The S=3/2 site
# is a custom site of dmrgpy, so all three have to give the same energy.
import numpy as np
from dmrgpy import spinchain

def energy(version, n):
    sc = spinchain.Spin_Chain(["S=3/2" for i in range(n)], itensor_version=version)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = 40 # bond dimension
    sc.nsweeps = 10 # number of DMRG sweeps
    return sc.gs_energy()

ns = [3, 4, 5] # chain lengths
e2 = [energy(2, n) for n in ns]
e3 = [energy(3, n) for n in ns]
epy = [energy("python", n) for n in ns]
print("Max difference v2-v3 =", np.max(np.abs(np.array(e2) - e3)))
print("Max difference v3-python =", np.max(np.abs(np.array(e3) - epy)))

import matplotlib.pyplot as plt
plt.plot(ns, e2, "o-", linewidth=4, label="ITensor v2")
plt.plot(ns, e3, "s--", linewidth=2, label="ITensor v3")
plt.plot(ns, epy, "x:", linewidth=2, label="python")
plt.xlabel("Chain length")
plt.ylabel("Ground state energy")
plt.legend()
plt.show()
