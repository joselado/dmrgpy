# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import fermionchain

def get(U=1.,mode="ED"):
    n = 8 # number of spinless fermionic sites
    fc = fermionchain.Fermionic_Chain(n) # create the chain
    h = 0
    for i in range(n-1): # hopping
        h = h + fc.Cdag[i]*fc.C[i+1]
        h = h + U*(fc.N[i]-.5)*(fc.N[i+1]-.5)
    for i in range(n-2): # hopping
        h = h + 0.1*fc.Cdag[i]*fc.C[i+2] # some small NNN
    h = h + 0.5*fc.N[0] # a local perturbation
    h = h + h.get_dagger()
    ##############################
    # Setup the Many Body Hamiltonian
    fc.set_hamiltonian(h) # set the hoppings
    wf = fc.get_gs(mode=mode) # get the ground state
    return wf.get_correlation_entropy_density() # compute the corr ent density

Sed = get(mode="ED") # get from ED
Stn = get(mode="DMRG") # get from MPS
inds = range(len(Sed)) # site indexes
import matplotlib.pyplot as plt

plt.plot(inds,Sed,label="ED",c="red")
plt.scatter(inds,Stn,label="DMRG",c="blue")
plt.legend()

plt.xlabel("Site")
plt.ylabel("Correlation entropy density")
plt.show()











