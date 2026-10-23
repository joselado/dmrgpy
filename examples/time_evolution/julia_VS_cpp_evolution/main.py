# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Cost of the ground state calculation of the S=1/2 Heisenberg chain
# as a function of the length, with the C++ ITensor backend and with the
# live Julia backend (ITensors.jl). The time is the DMRG wall time.
import time
from dmrgpy import spinchain

def get_time(n, julia):
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
    if julia: sc.setup_julia() # switch to the Julia backend
    h = 0
    for i in range(n-1): # Heisenberg Hamiltonian
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = 20 # bond dimension
    sc.nsweeps = 6 # number of DMRG sweeps
    t0 = time.time()
    sc.gs_energy()
    return time.time()-t0

ns = [10, 20, 40, 80] # chain lengths
get_time(10, True) # warm-up, Julia compiles the first call
tc = [get_time(n, False) for n in ns]
tj = [get_time(n, True) for n in ns]

import matplotlib.pyplot as plt
plt.plot(ns, tc, "o-", label="C++")
plt.plot(ns, tj, "s--", label="Julia")
plt.xlabel("number of sites")
plt.ylabel("time (s)")
plt.legend()
plt.show()
