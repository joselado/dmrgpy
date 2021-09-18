# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
from dmrgpy import spinfermionchain
n = 7
spins = [2 for i in range(n)] # spin 1/2 heisenberg chain
def geth(sc):
  h = 0
  for i in range(n-1): 
      h = h + sc.Sx[i]*sc.Sx[i+1]
      h = h + sc.Sy[i]*sc.Sy[i+1]
      h = h + sc.Sz[i]*sc.Sz[i+1]
  return h

sc = spinchain.Spin_Chain(spins) # create the spin chain
fc = spinfermionchain.Spin_Fermion_Hamiltonian(["S" for s in spins]) 

sc.set_hamiltonian(geth(sc)) # create Hamiltonian
fc.set_hamiltonian(geth(fc)) # create Hamiltonian

print("Energy with spins",sc.gs_energy())
print("Energy with fermions",fc.gs_energy())
print(fc.get_density())











