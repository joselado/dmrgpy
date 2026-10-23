# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Ground state energy per site of the S=1/2 Heisenberg chain, E_0/n,
# as a function of the chain length n. For an open chain it approaches
# the Bethe ansatz value of the infinite chain, 1/4 - ln(2).
import numpy as np
from dmrgpy import spinchain

ns = [4, 8, 12, 16, 20, 24, 30] # chain lengths
es = []
for n in ns:
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
    # Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = 40 # bond dimension
    sc.nsweeps = 10 # number of DMRG sweeps
    es.append(sc.gs_energy())
    print("n =", n, " E0 =", es[-1])

import matplotlib.pyplot as plt
plt.plot(ns, np.array(es)/np.array(ns), "o-", label="DMRG, open chain")
plt.axhline(0.25-np.log(2), color="black", linestyle="--", label="Bethe ansatz")
plt.xlabel("chain length n")
plt.ylabel("energy per site E0/n")
plt.legend()
plt.show()
