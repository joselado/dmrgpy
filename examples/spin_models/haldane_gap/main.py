# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Haldane gap of the S=1 Heisenberg chain. An open chain has four nearly
# degenerate lowest states, from the S=1/2 spins at the two edges. The bulk
# gap is then E_5 - E_1, and it tends to about 0.41 for a long chain.
import numpy as np
from dmrgpy import spinchain

gaps = []
ns = [8, 12, 16] # chain lengths
for n in ns:
    sc = spinchain.Spin_Chain(["S=1" for i in range(n)])
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = 20 # bond dimension
    sc.nsweeps = 20 # number of DMRG sweeps
    es = sc.get_excited(n=5) # five lowest states
    gaps.append(es[4] - es[0])
    print("n =", n, " gap =", gaps[-1])

import matplotlib.pyplot as plt
plt.plot(ns, gaps, "o-")
plt.axhline(0.41, color="gray", linestyle="--")
plt.xlabel("Chain length")
plt.ylabel("Haldane gap")
plt.show()
