# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
n = 5
fc = fermionchain.Fermionic_Hamiltonian(n) # create the chain
#m = np.matrix(np.random.random((2*n,2*n))) # random matrix
#m = m +m.T
def ft(i,j):
#    return m[i,j]
    if abs(j-i)==1: return 1.0 #+ np.random.random()
    return 0.0
fc.set_hoppings(ft) # hoppings
e0 = fc.gs_energy(mode="ED") # energy with exact diagonalization
e1 = fc.gs_energy(mode="DMRG") # energy with DMRG
print(e0,e1)
