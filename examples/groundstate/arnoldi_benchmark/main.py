# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Ground state with the Arnoldi algorithm, which only needs the action of
# the Hamiltonian on a state. We use it on a Heisenberg chain with random
# couplings, with the Hamiltonian applied as an MPO on matrix product
# states, and compare the energy with exact diagonalization (ED).
import numpy as np
from dmrgpy import spinchain
from dmrgpy.algebra import arnolditk

arnolditk.arnoldimode = "DMRG" # apply H on matrix product states
rng = np.random.RandomState(1)
ns = [4, 6, 8, 10] # chain lengths
eed, earn = [], []
for n in ns:
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
    h = 0
    for i in range(n-1): # random exchange couplings
        h = h + rng.random()*sc.Sx[i]*sc.Sx[i+1]
        h = h + rng.random()*sc.Sy[i]*sc.Sy[i+1]
        h = h + rng.random()*sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = 40 # bond dimension
    eed.append(sc.gs_energy(mode="ED"))
    sc.mode = None
    (e, wf) = arnolditk.mpsarnoldi(sc, h, mode="GS")
    earn.append(np.real(np.atleast_1d(e)[0]))
print("Max difference Arnoldi-ED =", np.max(np.abs(np.array(eed)-earn)))

import matplotlib.pyplot as plt
plt.plot(ns, eed, linewidth=4, label="ED")
plt.plot(ns, earn, "--", linewidth=2, label="Arnoldi")
plt.xlabel("number of sites")
plt.ylabel("ground state energy")
plt.legend()
plt.show()
