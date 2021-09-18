# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
from dmrgpy import multioperator
n = 8
fc = fermionchain.Fermionic_Chain(n,spinful=False) # create the chain
m = np.matrix(np.random.random((n,n)) + 1j*np.random.random((n,n)))
m = m + m.H
def ft(i,j):
    return m[i,j]
    if abs(j-i)==1: return 1.0 
    return 0.0
fc.set_hoppings(ft) # hoppings

den = 0

for i in range(n): # loop over sites
  deni = multioperator.obj2MO([["N",i]])
  den = den + deni

print("Total density with DMRG",fc.excited_vev(den,mode="DMRG",n=4).real)
print("Total density with ED",fc.excited_vev(den,mode="ED",n=4).real)

print("Energies with DMRG",fc.get_excited(mode="DMRG",n=4).real)
print("Energies with ED",fc.get_excited(mode="ED",n=4).real)









