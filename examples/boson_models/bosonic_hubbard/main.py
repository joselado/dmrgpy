# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Density profile <N_i> of three bosonic sites with onsite repulsion U,
# H = sum_i (b_i^dagger b_{i+1} + h.c.) + U sum_i n_i^2.
# A larger U makes extra bosons costly, so the density drops. This is
# what cold-atom experiments measure with in situ imaging.
# We use exact diagonalization.
import numpy as np
from dmrgpy import bosonchain

n = 3 # number of sites
Us = [0.0, 0.1, 0.2, 0.4] # onsite interaction

import matplotlib.pyplot as plt
for U in Us:
    bc = bosonchain.Bosonic_Chain(n)
    h = 0
    for i in range(n-1): # hopping
        h = h + bc.Adag[i]*bc.A[i+1] + bc.Adag[i+1]*bc.A[i]
    for i in range(n): # onsite repulsion
        h = h + U*bc.Adag[i]*bc.A[i]*bc.Adag[i]*bc.A[i]
    bc.set_hamiltonian(h)
    ds = [bc.vev(bc.N[i], mode="ED").real for i in range(n)]
    print("U =", U, " energy =", bc.gs_energy(mode="ED"))
    plt.plot(range(n), ds, "o-", label="U=%s" % U)
plt.xlabel("site i")
plt.ylabel("<N_i>")
plt.legend()
plt.show()
