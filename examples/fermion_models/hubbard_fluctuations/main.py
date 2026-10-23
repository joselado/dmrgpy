# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Local density fluctuation <N_i^2> - <N_i>^2 of the Hubbard chain at
# half filling, H = sum t c^dag c + U sum (n_up - 1/2)(n_dn - 1/2).
# For large U the electrons localize and the fluctuation goes to zero,
# and it is large for U=0. It is related to the double occupation.
import numpy as np
from dmrgpy import fermionchain

n = 8 # number of spinful sites
us = np.linspace(-10.0, 10.0, 11) # Hubbard interactions
ds = []
for u in us:
    fc = fermionchain.Spinful_Fermionic_Chain(n)
    h = 0
    for i in range(n-1): # hopping, spin up and spin down
        h = h + fc.Cdagup[i]*fc.Cup[i+1] + fc.Cdagdn[i]*fc.Cdn[i+1]
    h = h + h.get_dagger()
    for i in range(n): # Hubbard interaction
        h = h + u*(fc.Nup[i]-0.5)*(fc.Ndn[i]-0.5)
    fc.set_hamiltonian(h)
    fc.maxm = 20 # bond dimension
    fc.nsweeps = 6 # number of DMRG sweeps
    ds.append(np.mean(fc.get_density_fluctuation()))

import matplotlib.pyplot as plt
plt.plot(us, ds, marker="o")
plt.xlabel("Hubbard interaction U")
plt.ylabel("density fluctuation")
plt.show()
