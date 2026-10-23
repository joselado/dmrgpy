# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Fluctuation <Mz^2> of the total magnetization Mz = sum_i Sz_i in the
# ground state of the transverse-field Ising chain. We compute it with
# vev(Mz, npow=2), which applies Mz twice to the state, and by brute force
# with vev(Mz*Mz). Both agree, but the first one is faster for long chains.
import time
import numpy as np
from dmrgpy import spinchain

def get_times(n):
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
    h = 0
    for i in range(n-1): h = h - sc.Sz[i]*sc.Sz[i+1] # exchange
    for i in range(n): h = h + 0.6*sc.Sx[i] # transverse field
    sc.set_hamiltonian(h)
    sc.maxm = 40 # bond dimension
    sc.nsweeps = 10 # number of DMRG sweeps
    sc.get_gs()
    Mz = sum(sc.Sz)
    t0 = time.time()
    d1 = sc.vev(Mz, npow=2).real
    t1 = time.time()
    d2 = sc.vev(Mz*Mz).real
    t2 = time.time()
    print("n =", n, " difference =", abs(d1-d2))
    return t1-t0, t2-t1

ns = [4, 6, 8, 10, 12] # chain lengths
times = np.array([get_times(n) for n in ns])

import matplotlib.pyplot as plt
plt.plot(ns, times[:, 1]/times[:, 0], marker="o")
plt.xlabel("number of sites n")
plt.ylabel("time brute force / time npow=2")
plt.show()
