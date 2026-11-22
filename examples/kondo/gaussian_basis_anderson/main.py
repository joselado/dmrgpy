# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Anderson impurity on a uniform bath, solved with DMRG in two bases.
# In the site basis the Fermi sea of the bath is entangled everywhere.
# In the Gaussian basis the bath orbitals are rotated by the circuit of the
# U=0 ground state, so that the bath is almost a product state and DMRG
# needs a much smaller bond dimension (arXiv:2212.09798).
import numpy as np
from dmrgpy import fermionchain

n = 32 # orbitals: the impurity (orbital 0) and 31 bath sites
U = 4.0 # repulsion on the impurity
h = np.zeros((n,n)) # hopping matrix, the same for both spins
for i in range(n-1): h[i,i+1] = h[i+1,i] = -1.0

maxms = [8, 16, 32, 64] # bond dimensions
energies = {"site": [], "gaussian": []}
for basis in energies:
    for m in maxms:
        fc = fermionchain.Spinful_Fermionic_Chain(n)
        fc.maxm = m # bond dimension
        fc.nsweeps = 12 # number of DMRG sweeps
        H = 0
        for i in range(n):
            for j in range(n):
                if h[i,j] != 0: H = H + h[i,j]*(fc.Cdagup[i]*fc.Cup[j] + fc.Cdagdn[i]*fc.Cdn[j])
        H = H + U*(fc.Nup[0]-0.5)*(fc.Ndn[0]-0.5)
        if basis == "gaussian": # rotate the bath, the impurity stays as it is
            gb = fc.get_gaussian_basis(h=h, nf=n, fixed=[0])
            H = gb.rotate(H)
        fc.set_hamiltonian(H)
        energies[basis].append(fc.gs_energy())
    print(basis, "basis, energies:", np.round(energies[basis], 6))

# maxm=16 in the Gaussian basis is better than maxm=64 in the site basis
assert energies["gaussian"][1] < energies["site"][3]

import matplotlib.pyplot as plt
plt.plot(maxms, energies["site"], "o-", linewidth=3, label="site basis")
plt.plot(maxms, energies["gaussian"], "x--", markersize=10, label="Gaussian basis")
plt.xscale("log")
plt.xlabel("bond dimension maxm")
plt.ylabel("ground-state energy")
plt.legend()
plt.show()
