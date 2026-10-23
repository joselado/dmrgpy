# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Ground state of one chosen particle number N, with set_conserved_sector.
# Without it DMRG returns the global ground state, whatever N it has.
# E(N) of a spinless t-V chain, H = sum c^+_i c_i+1 + h.c. + V n_i n_i+1,
# gives the addition spectrum and the charge gap E(N+1)+E(N-1)-2E(N).
# We compute every sector with DMRG and compare with ED in the same sector.
import numpy as np
from dmrgpy import fermionchain

n = 10 # number of sites
V = 2.0 # nearest-neighbor repulsion
fc = fermionchain.Fermionic_Chain(n)
h = 0
for i in range(n-1):
    h = h + fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i] + V*fc.N[i]*fc.N[i+1]
fc.set_hamiltonian(h)
fc.maxm = 40 # bond dimension
fc.nsweeps = 10 # number of DMRG sweeps

Ns = list(range(n+1)) # particle numbers
e1, e2 = [], []
for N in Ns:
    fc.set_conserved_sector(Nf=N)
    e1.append(fc.gs_energy(mode="DMRG"))
    e2.append(fc.gs_energy(mode="ED"))
print("Max difference DMRG-ED =", np.max(np.abs(np.array(e1)-e2)))

import matplotlib.pyplot as plt
plt.plot(Ns, e2, "o", markersize=10, label="ED")
plt.plot(Ns, e1, "x", markersize=10, label="DMRG")
plt.xlabel("particle number N")
plt.ylabel("ground state energy E(N)")
plt.legend()
plt.show()
