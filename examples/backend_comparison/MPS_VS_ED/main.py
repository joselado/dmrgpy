# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain
n = 6
e0s,e1s = [],[]
for i in range(5):
#  spins = [np.random.randint(2,6) for i in range(n)] # spin 1/2 heisenberg chain
  spins = [np.random.randint(2,5) for i in range(n)] # spin 1/2 heisenberg chain
  sc = spinchain.Spin_Chain(spins) # create the spin chain
  ms = [np.random.random((3,3)) for i in range(n)]
  ms = [m + m.T for m in ms]
  h = 0
  Si = [sc.Sx,sc.Sy,sc.Sz]
  for i in range(n):
      for j in range(n):
          h = h + sc.Sx[i]*sc.Sx[j]*np.random.random()
          h = h + sc.Sy[i]*sc.Sy[j]*np.random.random()
          h = h + sc.Sz[i]*sc.Sz[j]*np.random.random()
  sc.hamiltonian = h
  e0 = sc.gs_energy(mode="DMRG") # compute the ground state energy
  e1 = sc.gs_energy(mode="ED") # compute the ground state energy
  print("Spin chain",spins)
  print("Energy with ED",e1)
  print("Energy with DMRG",e0)
  print("\n")
  de = np.abs(e1-e0)
  if de>0.1: raise
  e0s.append(e0)
  e1s.append(e1)
print("Test passed")

import matplotlib.pyplot as plt
lims = [min(e0s+e1s),max(e0s+e1s)]
plt.plot(lims,lims,color="gray",linestyle="--",label="ED = DMRG")
plt.scatter(e1s,e0s,marker="o")
plt.xlabel("Energy with ED")
plt.ylabel("Energy with DMRG")
plt.legend()
plt.show()











