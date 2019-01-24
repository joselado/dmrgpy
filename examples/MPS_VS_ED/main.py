# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
n = 4
for i in range(10):
  spins = [np.random.randint(2,6) for i in range(n)] # spin 1/2 heisenberg chain
  sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain
  def fj(i,j): 
      return 0.3
#      m = np.random.random((3,3)) # random couplings
#      return m + m.transpose()
  def fb(i): 
      m = np.random.random(3)
      return m + m.transpose()
  sc.set_exchange(fj)
  sc.set_fields(fb)
  e0 = sc.gs_energy(mode="DMRG") # compute the ground state energy
  e1 = sc.gs_energy(mode="ED") # compute the ground state energy
  print("Spin chain",spins)
  print("Energy with ED",e1)
  print("Energy with DMRG",e0)
  print("\n")
  de = np.abs(e1-e0)
  if de>0.1: raise
print("Test passed")


