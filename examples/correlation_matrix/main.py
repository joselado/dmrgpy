# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import fermionchain
n = 10 # number of spinful fermionic sites
fc = fermionchain.Spinful_Fermionic_Chain(n) # create the chain
h = 0
U = 0.0
for i in range(n-1): # hopping
    h = h + fc.Cdagup[i]*fc.Cup[i+1]
    h = h + fc.Cdagdn[i]*fc.Cdn[i+1]
for i in range(n): # Hubbard
    h = h + U*(fc.Nup[i]-.5)*(fc.Ndn[i]-.5)
h = h + h.get_dagger()
##############################
# Setup the Many Body Hamiltonian
fc.maxm = 40
fc.nsweeps = 10
fc.set_hamiltonian(h) # set the hoppings
fc.get_gs()

cm = fc.get_correlation_matrix()
import scipy.linalg as lg
print(lg.eigvalsh(cm))
