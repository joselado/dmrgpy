# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Correlated orbitals of a spinful chain with a Hubbard interaction U on
# the first site only. Diagonalizing the correlation matrix <c^dag_i c_j>
# gives orbitals that are either full, empty or entangled with the rest.
# out[i][j] is the weight of orbital i on site j.
import numpy as np
from dmrgpy import fermionchain

n = 12 # number of spinful sites
U = 2.0 # Hubbard interaction on site 0
fc = fermionchain.Spinful_Fermionic_Chain(n)
h = 0
for i in range(n-1): # hopping
    h = h + fc.Cdagup[i]*fc.Cup[i+1] + fc.Cdagdn[i]*fc.Cdn[i+1]
h = h + U*(fc.Nup[0]-.5)*(fc.Ndn[0]-.5) # Hubbard on site 0
h = h + h.get_dagger()
fc.maxm = 40 # bond dimension
fc.nsweeps = 10 # number of DMRG sweeps
fc.set_hamiltonian(h)
fc.get_gs()
out = np.array(fc.get_correlated_orbitals(ordered=False))

import matplotlib.pyplot as plt
plt.imshow(out, aspect="auto", origin="lower")
plt.colorbar(label="orbital weight")
plt.xlabel("site j")
plt.ylabel("orbital i")
plt.show()
