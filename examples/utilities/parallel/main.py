# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Ground-state energy of a Heisenberg chain with a staggered field u,
# H = sum_i S_i . S_{i+1} + u sum_i (-1)^i Sz_i, for several values of u.
# Each value of u is independent, so we solve them in parallel with
# parallel.pcall, which runs the function on several cores.
import numpy as np
from dmrgpy import spinchain, parallel

def get(u):
    n = 16 # number of sites
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(n):
        h = h + u*(-1)**i*sc.Sz[i] # staggered field
    sc.set_hamiltonian(h)
    sc.maxm = 40 # bond dimension
    sc.nsweeps = 10 # number of DMRG sweeps
    return sc.gs_energy()

parallel.cores = 4 # number of cores
us = np.linspace(0, 1, 4) # values of the staggered field
es = parallel.pcall(get, us)

import matplotlib.pyplot as plt
plt.plot(us, es, marker="o")
plt.xlabel("staggered field u")
plt.ylabel("ground-state energy")
plt.show()
