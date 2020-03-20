# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
from dmrgpy import multioperator
n = 6
fc = fermionchain.Fermionic_Chain(n,spinful=False) # create the chain
m = np.matrix(np.random.random((n,n)) + 1j*np.random.random((n,n)))
m = m + m.H
def ft(i,j):
    return m[i,j]
    if abs(j-i)==1: return 1.0 
    return 0.0

h = 0
for i in range(n):
  for j in range(n):
    h = h + fc.Cdag[i]*fc.C[j]*ft(i,j)


fc.set_hamiltonian(h) # hoppings

print("Energy with ED",fc.gs_energy(mode="ED"))
print("Energy with DMRG",fc.gs_energy(mode="DMRG"))

pairs = [(0,i) for i in range(n)]

den = []

for i in range(n): # loop over sites
  cd = multioperator.obj2MO([["Cdag",i]])
  c = multioperator.obj2MO([["C",i]])
  den.append(cd*c)

print("Density DMRG",[fc.vev(di,mode="DMRG").real for di in den])
print("Density ED",[fc.vev(di,mode="ED").real for di in den])

