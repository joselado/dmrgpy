# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Addition spectrum of a spinless t-V chain. We compute the ground state
# energy E(N) inside each sector of fixed particle number N, with exact
# diagonalization. The charge gap, Delta(N) = E(N+1) + E(N-1) - 2 E(N),
# is the energy to add and remove a particle, and it is measured in
# transport and tunneling experiments.
import numpy as np
from dmrgpy import fermionchain

n = 8 # number of sites
V = 2.0 # nearest neighbor repulsion
fc = fermionchain.Fermionic_Chain(n)
h = 0
for i in range(n-1):
    h = h + fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i]
    h = h + V*fc.N[i]*fc.N[i+1]
fc.set_hamiltonian(h)
e_all = fc.gs_energy(mode="ED") # ground state without any constraint

nfs = np.arange(n+1) # particle numbers
es = []
for nf in nfs:
    fc.set_conserved_sector(Nf=int(nf)) # only states with nf particles
    es.append(fc.gs_energy(mode="ED"))
es = np.array(es)
assert abs(np.min(es) - e_all) < 1e-8 # the lowest sector is the ground state
gaps = es[2:] + es[:-2] - 2*es[1:-1] # charge gap

import matplotlib.pyplot as plt
plt.plot(nfs[1:-1], gaps, "o-")
plt.xlabel("Number of particles N")
plt.ylabel("E(N+1) + E(N-1) - 2 E(N)")
plt.show()
