# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
n = 16
fc = fermionchain.Fermionic_Chain(n) # create the chain
def ft(i,j):
    if abs(j-i)==1: return 1.0 #+ np.random.random()
    return 0.0
fc.set_hoppings(ft) # hoppings
e0 = fc.gs_energy(mode="ED") # energy with exact diagonalization
e1 = fc.gs_energy(mode="DMRG") # energy with DMRG
print("Energy with ED",e0)
print("Energy with DMRG",e1)


