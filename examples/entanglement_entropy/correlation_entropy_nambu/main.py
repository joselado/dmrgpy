# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import fermionchain

def get(U):
    n = 20 # number of spinful fermionic sites
    fc = fermionchain.Fermionic_Chain(n) # create the chain
    h = 0
    for i in range(n-1): # hopping
        h = h + fc.Cdag[i]*fc.C[i+1]
        h = h + U*(fc.N[i]-.5)*(fc.N[i+1]-.5)
#    for i in range(n-1): # superconductivity
#        h = h + .4*fc.C[i]*fc.C[i+1]
    h = h + h.get_dagger()
    ##############################
    # Setup the Many Body Hamiltonian
    fc.set_hamiltonian(h) # set the hoppings
    wf = fc.get_gs(mode="DMRG")
    return wf.get_correlation_entropy(basis="Nambu")


print(get(0.2))
exit()

Us = np.linspace(-2.,0.,60)
Ss = [get(U) for U in Us]

import matplotlib.pyplot as plt

np.savetxt("ENTROPY.OUT",np.array([Us,Ss]).T)
plt.plot(Us,Ss)
plt.show()











