# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
n = 10
fc = fermionchain.Fermionic_Hamiltonian(n,spinful=False) # create the chain
def ft(i,j):
    if abs(j-i)==1: return 1.0 
    if i==j: return 0.5
    return 0.0
fc.set_hoppings(ft) # hoppings

e = fc.gs_energy() # energy with DMRG
print("Ground state energy",e)


import matplotlib.pyplot as plt

# Now compute the ground state correlator
pairs = [(0,i) for i in range(n)]
cs0 = fc.get_correlator(pairs=pairs,mode="DMRG")
cs1 = fc.get_correlator(pairs=pairs,mode="ED")
x = np.array(range(len(cs0))) # distances
plt.plot(x,cs0,c="blue",label="DMRG")
plt.scatter(x,cs1,c="red",label="ED")
plt.legend()
plt.show()


