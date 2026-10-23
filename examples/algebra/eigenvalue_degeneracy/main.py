# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Ground state degeneracy of the ferromagnetic Heisenberg dimer.
# The ground state is the S=1 triplet, so the degeneracy is 3.
# It is estimated from the weight of the pole of a correlator at E0,
# and it should not depend on the broadening delta of the poles.
import numpy as np
from dmrgpy import spinchain
from dmrgpy.degeneracy import pole_eigenvalue_degeneracy

n = 2 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])

# ferromagnetic Heisenberg Hamiltonian, H = -sum_i S_i . S_{i+1}
h = 0
for i in range(n-1):
    h = h - sc.Sx[i]*sc.Sx[i+1] - sc.Sy[i]*sc.Sy[i+1] - sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 20 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps
e0 = sc.gs_energy()
print("Ground state energy =", e0)

deltas = [0.05, 0.1, 0.2, 0.3, 0.5] # broadenings of the poles
degs = [float(pole_eigenvalue_degeneracy(sc, h, e0, delta=d).real) for d in deltas]
print("Degeneracies =", degs)

import matplotlib.pyplot as plt
plt.plot(deltas, degs, "o-", label="estimated")
plt.axhline(3, color="gray", linestyle="--", label="exact")
plt.xlabel("broadening delta")
plt.ylabel("ground state degeneracy")
plt.legend()
plt.show()
