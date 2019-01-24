# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
n = 40
fc = fermionchain.Fermionic_Hamiltonian(n,spinful=False) # create the chain
def ft(i,j):
    if abs(j-i)==1: return 1.0 
    return 0.0
fc.set_hoppings(ft) # hoppings
fc.set_hubbard(lambda i,j: ft(i,j)*0.) # density-density interaction

e = fc.gs_energy() # energy with DMRG
print("Ground state energy",e)


import matplotlib.pyplot as plt

# Now compute the ground state correlator
pairs = [(0,i) for i in range(n-1)]
cs = fc.get_correlator(pairs=pairs)
x = np.array(range(len(cs))) # distances
plt.plot(x,cs,c="blue",label="DMRG")
plt.xlabel("N")
plt.xlabel("Correlator")
plt.legend()
plt.show()


