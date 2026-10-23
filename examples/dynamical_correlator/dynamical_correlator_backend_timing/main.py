# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Cost of the two main ways to get a local spin spectrum S(w): the kernel
# polynomial method (KPM, Chebyshev moments) and the time-dependent method
# (TD, real-time evolution and a Fourier transform). We time both as a
# function of the chain length of the S=1/2 Heisenberg chain.
import time
import numpy as np
from dmrgpy import spinchain

es = np.linspace(-1.0, 6.0, 300) # energies
delta = 0.3 # broadening of the peaks
ns = [8, 12, 16] # chain lengths
times = {"KPM": [], "TD": []}
for n in ns:
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = 40 # bond dimension
    sc.nsweeps = 10 # number of DMRG sweeps
    sc.gs_energy() # ground state, shared by both methods
    name = (sc.Sz[n//2], sc.Sz[n//2])
    for submode in times:
        t0 = time.time()
        sc.get_dynamical_correlator(name=name, es=es, delta=delta, submode=submode)
        times[submode].append(time.time() - t0)
    print("n =", n, "times (s) =", times["KPM"][-1], times["TD"][-1])

import matplotlib.pyplot as plt
for submode in times:
    plt.plot(ns, times[submode], "o-", label=submode)
plt.xlabel("chain length n")
plt.ylabel("wall time (s)")
plt.legend()
plt.show()
