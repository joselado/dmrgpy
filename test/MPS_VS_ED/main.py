import sys
import os
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg
import spinchain

n = 5
for i in range(10):
  spins = [np.random.randint(2,7) for i in range(n)] # spin 1/2 heisenberg chain
  sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain
  e0 = sc.gs_energy() # compute the ground state energy
#  print("Energy per site, DMRG",e0/n)
  e1 = sc.gs_energy(mode="full") # compute the ground state energy
  de = np.abs(e1-e0)
  if de>0.0001: raise

print("Test passed")
