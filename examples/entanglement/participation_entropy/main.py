# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Second Renyi participation entropy S2 = -log sum_s p(s)^2 of the
# Heisenberg ground state in the Sz basis, with p(s) = |psi(s)|^2 the
# probability of a snapshot s. DMRG builds p as an MPS by recursive
# sketched interpolation (arXiv:2602.17974). S2 grows linearly with the
# length; we check it against exact diagonalization (ED).
import numpy as np
from dmrgpy import spinchain

ns = [4, 6, 8, 10, 12] # chain lengths
s2, s2ed = [], []
for n in ns:
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
    sc.maxm = 40 # bond dimension
    sc.nsweeps = 10 # number of DMRG sweeps
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    s2.append(sc.get_participation_entropy(q=2, maxm=60))
    s2ed.append(sc.get_participation_entropy(q=2, mode="ED"))
print("Max difference DMRG-ED =", np.max(np.abs(np.array(s2)-s2ed)))
assert np.max(np.abs(np.array(s2)-s2ed)) < 1e-4

import matplotlib.pyplot as plt
plt.plot(ns, s2ed, linewidth=4, label="ED")
plt.plot(ns, s2, "o--", label="DMRG + RSI")
plt.xlabel("number of sites")
plt.ylabel("S2 (participation)")
plt.legend()
plt.show()
