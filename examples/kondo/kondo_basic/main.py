# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Kondo impurity: a spin on site 0 (with a Hubbard term, so that it stays
# half filled) coupled by the exchange JK to the end of a chain of
# spinful fermions. We compute the correlation entropy of every orbital,
# which measures how entangled each orbital is with the rest. The
# entanglement of the impurity and of the first chain site grows with JK.
import numpy as np
from dmrgpy import fermionchain

def entropy(JK):
    n = 8 # number of sites
    fc = fermionchain.Spinful_Fermionic_Chain(n)
    h = 0
    for i in range(1, n-1): # hopping along the chain
        h = h + fc.Cdagup[i]*fc.Cup[i+1] + fc.Cdagdn[i]*fc.Cdn[i+1]
    h = h + (fc.Nup[0]-.5)*(fc.Ndn[0]-.5) # Hubbard term in the impurity
    h = h + JK/2.*(fc.Sx[0]*fc.Sx[1] + fc.Sy[0]*fc.Sy[1] + fc.Sz[0]*fc.Sz[1])
    h = h + h.get_dagger() # make it Hermitian
    fc.set_hamiltonian(h)
    fc.maxm = 40 # bond dimension
    fc.nsweeps = 10 # number of DMRG sweeps
    return fc.get_gs().get_correlation_entropy_density()

import matplotlib.pyplot as plt
for JK in [0.0, 0.5, 1.0, 2.0]:
    s = entropy(JK)
    plt.plot(range(len(s)), s, "o-", label="JK=%s" % JK)
plt.xlabel("orbital")
plt.ylabel("correlation entropy")
plt.legend()
plt.show()
