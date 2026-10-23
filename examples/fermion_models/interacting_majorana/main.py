# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Local density of states of a chain of spinless fermions with
# nearest neighbor interaction V, H = sum_i (c_i^dagger c_{i+1} + h.c.)
# + V sum_i (n_i - 1/2)(n_{i+1} - 1/2),
# DOS(w) = sum_n |<n|c^dagger_0|GS>|^2 delta(w - E_n + E_0).
# It is what STM measures at the end of the chain. We use DMRG with KPM.
import numpy as np
from dmrgpy import fermionchain

n = 20 # number of sites
V = -0.9 # interaction
fc = fermionchain.Fermionic_Chain(n)
h = 0
for i in range(n-1):
    h = h + fc.Cdag[i]*fc.C[i+1]
    h = h + V*(fc.N[i]-0.5)*(fc.N[i+1]-0.5)
h = h + h.get_dagger()
fc.set_hamiltonian(h)
fc.maxm = 20 # bond dimension
fc.nsweeps = 10 # number of DMRG sweeps

i = 0 # site at the edge
es = np.linspace(0., 3., 300) # energies
(x, y) = fc.get_dynamical_correlator(es=es, name=(fc.Cdag[i], fc.C[i]), delta=0.1)

import matplotlib.pyplot as plt
plt.plot(x, y.real)
plt.xlabel("energy w")
plt.ylabel("DOS at site 0")
plt.show()
