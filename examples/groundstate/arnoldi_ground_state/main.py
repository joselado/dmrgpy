# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Ground state energy of Heisenberg chains with random exchange, from
# an Arnoldi iteration on matrix product states, compared with DMRG.
# Both methods minimize the energy <GS|H|GS>, so they must agree.
import numpy as np
from dmrgpy import spinchain
from dmrgpy import mpsalgebra

np.random.seed(1)
ns = [4, 6, 8, 10] # chain lengths
e_arnoldi, e_dmrg = [], []
for n in ns:
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
    h = 0
    for i in range(n-1): # random first-neighbor exchange
        h = h + np.random.random()*sc.Sx[i]*sc.Sx[i+1]
        h = h + np.random.random()*sc.Sy[i]*sc.Sy[i+1]
        h = h + np.random.random()*sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    e_arnoldi.append(np.real(mpsalgebra.lowest_energy_arnoldi(sc, h)[0]).ravel()[0])
    e_dmrg.append(sc.gs_energy())
print("Max difference Arnoldi-DMRG =", np.max(np.abs(np.array(e_arnoldi)-e_dmrg)))

import matplotlib.pyplot as plt
plt.plot(ns, e_dmrg, "o-", linewidth=4, label="DMRG")
plt.plot(ns, e_arnoldi, "x--", linewidth=2, label="Arnoldi")
plt.xlabel("number of sites")
plt.ylabel("ground state energy")
plt.legend()
plt.show()
