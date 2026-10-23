# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Ground-state energy of the S=1/2 Heisenberg chain, H = sum_i S_i . S_{i+1},
# written in two ways: with spin operators, and with spins built from
# a fermion per site (spin-fermion chain). Both describe the same
# physics, so the energies must agree for every chain length.
import numpy as np
from dmrgpy import spinchain, spinfermionchain

def energy(chain, n):
    h = 0
    for i in range(n-1):
        h = h + chain.Sx[i]*chain.Sx[i+1] + chain.Sy[i]*chain.Sy[i+1] + chain.Sz[i]*chain.Sz[i+1]
    chain.set_hamiltonian(h)
    return chain.gs_energy()

ns = [3, 4, 5, 6, 7] # chain lengths
e1 = [energy(spinchain.Spin_Chain(["S=1/2" for i in range(n)]), n) for n in ns]
e2 = [energy(spinfermionchain.Spin_Fermion_Hamiltonian(["S" for i in range(n)]), n) for n in ns]
print("Max difference =", np.max(np.abs(np.array(e1)-e2)))

import matplotlib.pyplot as plt
plt.plot(ns, e1, "o", markersize=10, label="spin chain")
plt.plot(ns, e2, "x", markersize=10, label="spin-fermion chain")
plt.xlabel("number of sites")
plt.ylabel("ground-state energy")
plt.legend()
plt.show()
