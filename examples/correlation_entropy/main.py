# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import fermionchain

def get(U):
    n = 10 # number of spinful fermionic sites
    fc = fermionchain.Fermionic_Chain(n) # create the chain
    h = 0
    for i in range(n-1): # hopping
        h = h + fc.Cdag[i]*fc.C[i+1]
        h = h + U*(fc.N[i]-.5)*(fc.N[i+1]-.5)
    h = h + h.get_dagger()
    ##############################
    # Setup the Many Body Hamiltonian
    fc.maxm = 40
    fc.nsweeps = 10
    fc.set_hamiltonian(h) # set the hoppings
    fc.get_gs()
    return fc.get_correlation_entropy()

Us = np.linspace(-2.,0.,60)
Ss = [get(U) for U in Us]

np.savetxt("ENTROPY.OUT",np.array([Us,Ss]).T)
