# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Magnetization Mz = sum_i <Sz_i> of the S=1/2 Heisenberg chain versus the
# applied field B, H = sum_i S_i.S_{i+1} + B sum_i Sz_i. Each field is an
# independent calculation, so we run them in parallel with pcall (only
# on a Slurm cluster such as Triton). In serial, use [getm(b) for b in bs].
import numpy as np
from dmrgpy import spinchain, parallelslurm

def getm(b): # everything is defined inside, as pcall requires
    n = 100 # number of sites
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    Mz = 0
    for i in range(n):
        h = h + b*sc.Sz[i]
        Mz = Mz + sc.Sz[i]
    sc.set_hamiltonian(h)
    sc.maxm = 40 # bond dimension
    sc.nsweeps = 10 # number of DMRG sweeps
    return sc.vev(Mz, mode="DMRG").real

bs = np.linspace(0., 1., 40) # magnetic fields
ms = parallelslurm.pcall(getm, bs, time=1) # one task per field

import matplotlib.pyplot as plt
plt.plot(bs, ms, marker="o")
plt.xlabel("magnetic field B")
plt.ylabel("magnetization Mz")
plt.show()
