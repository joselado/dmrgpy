# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Ground-state energy of the S=1/2 Heisenberg chain with DMRG, computed
# with the C++ backend (ITensor v2) and with the Julia backend. Both must
# give the same energy, and we compare them as a function of chain length.
# Needs Julia and juliacall installed.
import numpy as np
from dmrgpy import spinchain

def get_energy(n, julia):
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
    if julia: sc.setup_julia() # Julia backend
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = 40 # bond dimension
    sc.nsweeps = 10 # number of DMRG sweeps
    return sc.gs_energy()

ns = [10, 14, 18, 22]
get_energy(8, True) # first call, Julia compiles
e_cpp = np.array([get_energy(n, False) for n in ns])
e_jul = np.array([get_energy(n, True) for n in ns])
print("Max difference =", np.max(np.abs(e_cpp - e_jul)))

import matplotlib.pyplot as plt
plt.plot(ns, e_cpp/np.array(ns), linewidth=4, label="C++")
plt.plot(ns, e_jul/np.array(ns), "--", linewidth=2, label="Julia")
plt.xlabel("chain length n")
plt.ylabel("energy per site")
plt.legend()
plt.show()
