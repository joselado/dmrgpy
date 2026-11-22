# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Anderson impurity on a uniform bath, solved with DMRG in two Gaussian bases.
# The first one is the circuit of the ground state without U.
# The second one is the circuit of the natural orbitals of the interacting
# ground state, the eigenvectors of its correlation matrix <c_i^dag c_j>.
# The second basis needs a smaller bond dimension for the same energy.
import numpy as np
from dmrgpy import fermionchain

n = 32 # orbitals: the impurity (orbital 0) and 31 bath sites
U = 4.0 # repulsion on the impurity
h = np.zeros((n,n)) # hopping matrix, the same for both spins
for i in range(n-1): h[i,i+1] = h[i+1,i] = -1.0
fc = fermionchain.Spinful_Fermionic_Chain(n)
fc.nsweeps = 12 # number of DMRG sweeps
H = 0
for i in range(n):
    for j in range(n):
        if h[i,j] != 0: H = H + h[i,j]*(fc.Cdagup[i]*fc.Cup[j] + fc.Cdagdn[i]*fc.Cdn[j])
H = H + U*(fc.Nup[0]-0.5)*(fc.Ndn[0]-0.5)

bases = {"U=0 state": fc.get_gaussian_basis(h=h, nf=n, fixed=[0])}
fc.set_gaussian_basis(bases["U=0 state"], H)
fc.maxm = 32 # the state the natural orbitals come from
fc.gs_energy()
bases["natural orbitals"] = fc.gaussian_basis.natural_basis()

maxms = [8, 12, 16, 24] # bond dimensions
energies = {name: [] for name in bases}
for name in bases:
    for m in maxms:
        fc.set_gaussian_basis(bases[name], H) # cold start in this basis
        fc.maxm = m
        energies[name].append(fc.gs_energy())
    print(name, "energies:", np.round(energies[name], 6))

# at maxm=8 the natural orbitals give a lower energy
assert energies["natural orbitals"][0] < energies["U=0 state"][0]

import matplotlib.pyplot as plt
plt.plot(maxms, energies["U=0 state"], "o-", linewidth=3, label="circuit of the U=0 state")
plt.plot(maxms, energies["natural orbitals"], "x--", markersize=10, label="natural orbitals")
plt.xlabel("bond dimension maxm")
plt.ylabel("ground-state energy")
plt.legend()
plt.show()
