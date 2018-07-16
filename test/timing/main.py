from __future__ import print_function
import sys
import os
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg
import spinchain

ns = [10,100,1000]
for n in ns:
#  spins = [np.random.randint(2,7) for i in range(n)] # spin 1/2 heisenberg chain
  spins = [2 for i in range(n)] # spin 1/2 heisenberg chain
  sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain
  e0 = sc.gs_energy(mode="DMRG") # compute the ground state energy
#  e1 = sc.gs_energy(mode="ED") # compute the ground state energy
#  print("Spin chain",spins)
#  print("Energy with ED",e1)
  print("Energy with DMRG",e0)
  print("\n")
#  de = np.abs(e1-e0)
#  if de>0.0001: raise

print("Test passed")
