# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Cost of a dynamical correlator computed with real-time evolution
# (submode "TD") as a function of the chain length n. The correlator is
# <Sz_0(t) Sz_1(0)> of a Heisenberg chain, Fourier transformed to energy.
# We fit the time to a power law, time ~ n^p, and print the exponent p.
import time
import numpy as np
from dmrgpy import spinchain

ns = [4, 6, 8, 10] # chain lengths
ts = []
for n in ns:
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = 10 # bond dimension
    sc.nsweeps = 2 # number of DMRG sweeps
    sc.get_gs()
    t0 = time.time()
    sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[1]), submode="TD",
            es=np.linspace(-1.0, 6.0, 400), dt=0.1, nt=200)
    ts.append(time.time() - t0)
p = np.polyfit(np.log(ns), np.log(ts), deg=1) # straight line in log-log
print("Exponent p =", p[0])

import matplotlib.pyplot as plt
plt.plot(np.log(ns), np.log(ts), "o", markersize=10, label="TD")
plt.plot(np.log(ns), np.polyval(p, np.log(ns)), "--", label="fit, p = %.2f" % p[0])
plt.xlabel("log n")
plt.ylabel("log time")
plt.legend()
plt.show()
