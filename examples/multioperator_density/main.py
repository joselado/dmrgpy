# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
from dmrgpy import multioperator
n = 4
fc = fermionchain.Fermionic_Hamiltonian(n,spinful=False) # create the chain
m = np.matrix(np.random.random((n,n)) + 1j*np.random.random((n,n)))
m = m + m.H
def ft(i,j):
    return m[i,j]
    if abs(j-i)==1: return 1.0 
    return 0.0
fc.set_hoppings(ft) # hoppings
pairs = [(0,i) for i in range(n)]

den = multioperator.zero()

for i in range(n): # loop over sites
  cd = multioperator.obj2MO([["Cdag",i]])
  c = multioperator.obj2MO([["C",i]])
  den = den + cd*c

print("Total density",fc.vev(den).real)

