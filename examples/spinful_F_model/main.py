# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import fermionchain

# example for a model containing spinful fermions
# and spinless pseudofermions in each site
n = 4 # number of fermionic sites (each having spin and f fermions)
fc = fermionchain.Spinful_F_Fermionic_Chain(n) # create the chain
# initialize Hamiltonian #
h = 0

# first make a Hubbard model for the spinful fermions
for i in range(n-1): # hopping
    h = h + fc.Cdagup[i]*fc.Cup[i+1]
    h = h + fc.Cdagdn[i]*fc.Cdn[i+1]
for i in range(n): # Hubbard
    h = h + (fc.Nup[i]-.5)*(fc.Ndn[i]-.5)

# now add kinetic energy for the pseudofermions
for i in range(n-1): # hopping
    h = h + fc.Fdag[i]*fc.F[i+1]


# and now add hybridization between spinful and pseudo fermions
for i in range(n): # local hybridization
    h = h + fc.Fdag[i]*(fc.Nup[i]+fc.Ndn[i]) # fermion times density

h = h + h.get_dagger() # add the Hermitian conjugate
##############################
# Setup the Many Body Hamiltonian
fc.set_hamiltonian(h) # set the hoppings
print(fc.gs_energy()) # compute the ground state energy
