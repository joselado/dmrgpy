# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Heisenberg chain of S=1 spins with S=1/2 spins at the two ends and a
# single-ion anisotropy, H = sum_i S_i.S_{i+1} - D sum_i (Sz_i)^2. We
# follow the entanglement entropy of the ground state as D grows. The
# DMRG starts from a Neel state, which is the ground state of a staggered
# field, so that it does not get stuck in the wrong sector.
import numpy as np
from dmrgpy import spinchain

n = 20 # number of S=1 spins
spins = ["S=1/2"] + ["S=1" for i in range(n)] + ["S=1/2"]
N = len(spins)

def entropy(D):
    sc = spinchain.Spin_Chain(spins)
    sc.nsweeps = 10 # number of DMRG sweeps
    sc.maxm = 20 # bond dimension
    h0 = 0 # staggered field, its ground state is the Neel state
    for i in range(N): h0 = h0 + (-1)**i*sc.Sz[i]
    sc.set_hamiltonian(h0)
    wf0 = sc.get_gs()
    h = 0
    for i in range(N-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(N): h = h - D*sc.Sz[i]*sc.Sz[i]
    sc.set_hamiltonian(h)
    wf = sc.get_gs(wf0=wf0) # start from the Neel state
    return wf.get_entropy(N//2) # entropy of the half chain

Ds = np.linspace(0.0, 0.6, 10)
ss = [entropy(D) for D in Ds]

import matplotlib.pyplot as plt
plt.plot(Ds, ss, marker="o")
plt.xlabel("single-ion anisotropy D")
plt.ylabel("entanglement entropy")
plt.show()
