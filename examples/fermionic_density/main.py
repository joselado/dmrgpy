# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
n = 4
fc = fermionchain.Fermionic_Hamiltonian(n,spinful=False) # create the chain
m = np.matrix(np.random.random((n,n)) + 1j*np.random.random((n,n)))
m = m + m.H
def ft(i,j):
    return m[i,j]
    if abs(j-i)==1: return 1.0 
    return 0.0
fc.set_hoppings(ft) # hoppings
e0 = fc.gs_energy(mode="ED") # energy with exact diagonalization
e1 = fc.gs_energy(mode="DMRG") # energy with DMRG
print("Energy with ED",e0)
print("Energy with DMRG",e1)

d1 = fc.get_density_fluctuation(mode="DMRG")
d2 = fc.get_density_fluctuation(mode="ED")

plt.scatter(range(len(d1)),d1,c="red",label="DMRG")
plt.plot(range(len(d2)),d2,c="blue",label="ED")
plt.legend()
plt.show()
