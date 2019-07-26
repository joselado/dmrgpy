# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
n = 8
fc = fermionchain.Fermionic_Hamiltonian(n,spinful=False) # create the chain
m = np.matrix(np.random.random((n,n)) + 1j*np.random.random((n,n)))
m = m + m.H
def ft(i,j):
#    return m[i,j]
    if abs(j-i)==1: return -1.0 
    return 0.0
fc.set_hoppings(ft) # hoppings
fc.maxm = 50
e0 = fc.gs_energy(mode="ED") # energy with exact diagonalization
e1 = fc.gs_energy(mode="DMRG") # energy with DMRG
print("Energy with ED",e0)
print("Energy with DMRG",e1)


