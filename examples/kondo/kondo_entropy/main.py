# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Kondo coupling of a spin to a short spinful fermion chain. Site 0 is a
# Hubbard site that holds the spin, and it is coupled with JK to site 1 of
# the chain. The correlation entropy of the ground state measures how far it
# is from a single Slater determinant, and it jumps when the Kondo coupling is switched on.
import numpy as np
from dmrgpy import fermionchain

def get_entropy(JK):
    n = 4 # number of sites
    fc = fermionchain.Spinful_Fermionic_Chain(n)
    h = 0
    for i in range(1, n-1): # hopping along the chain
        h = h + fc.Cdagup[i]*fc.Cup[i+1] + fc.Cdagdn[i]*fc.Cdn[i+1]
    h = h + (fc.Nup[0]-.5)*(fc.Ndn[0]-.5) # Hubbard, to keep half filling
    h = h + JK/2.*(fc.Sx[0]*fc.Sx[1] + fc.Sy[0]*fc.Sy[1] + fc.Sz[0]*fc.Sz[1])
    h = h + h.get_dagger()
    fc.maxm = 40 # bond dimension
    fc.set_hamiltonian(h)
    return fc.get_gs().get_correlation_entropy()

Js = np.linspace(0., 2., 8) # Kondo couplings
ss = [get_entropy(J) for J in Js]

import matplotlib.pyplot as plt
plt.plot(Js, ss, "o-")
plt.xlabel("Kondo coupling JK")
plt.ylabel("correlation entropy")
plt.show()
