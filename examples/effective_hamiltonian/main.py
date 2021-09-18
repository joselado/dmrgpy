# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import fermionchain
n = 2 # number of spinful fermionic sites
fc = fermionchain.Spinful_Fermionic_Chain(n) # create the chain
# initialize Hamiltonian #
h = 0
phi = 0.1
t = 1.0*np.exp(1j*phi*np.pi)
for i in range(n-1): # hopping
    h = h + t*fc.Cdagup[i]*fc.Cup[i+1]
    h = h + np.conjugate(t)*fc.Cdagdn[i]*fc.Cdn[i+1]
for i in range(n): # Hubbard
    h = h + 20*(fc.Nup[i]-.5)*(fc.Ndn[i]-.5)
h = h + h.get_dagger()

fc.set_hamiltonian(h)

# print the effective Hamiltonian in latex form
from dmrgpy import effectivehamiltonian
l = effectivehamiltonian.get_effective_hamiltonian(fc,method="single",
        mode="ED")
print("Effective Hamiltonian in latex form")
print(l) # write the Hamiltonian in latex












