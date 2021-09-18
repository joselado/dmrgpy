# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
n = 2
fc = fermionchain.Majorana_Chain(n) # create the chain

for i in range(n):
  for j in range(n):
    h = fc.G[i]*fc.G[j]*np.random.random()

h = h + h.get_dagger()

fc.set_hamiltonian(h)
e0 = fc.gs_energy(mode="ED") # energy with exact diagonalization

h = h.simplify()
fc.set_hamiltonian(h)
e1 = fc.gs_energy(mode="DMRG") # energy with DMRG
print("Energy with ED",e0)
print("Energy with DMRG",e1)











