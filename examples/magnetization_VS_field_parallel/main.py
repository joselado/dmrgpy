# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain

n = 10 # number of sites in your chain
spins = ["S=1/2" for i in range(n)] # create the sites
sc = spinchain.Spin_Chain(spins) # create the chain

# now define the Hamiltonian
def geth(b):
  h = 0
  for i in range(n-1): 
      h = h + sc.Sx[i]*sc.Sx[i+1] # add exchange
      h = h + sc.Sy[i]*sc.Sy[i+1] # add exchange
      h = h + sc.Sz[i]*sc.Sz[i+1] # add exchange
  for i in range(n): h = h + b*sc.Sz[i] # add transverse field
  return h

Mz = 0
for i in range(n): Mz = Mz + sc.Sz[i]

# function that you want to parallelize
def getm(b):
  h = geth(b)
  sc.set_hamiltonian(h) # and initialize the Hamiltonian
  mz = sc.vev(Mz,mode="DMRG").real
  return mz

bs = np.linspace(0.,1.,30) # values for which you want to call the function

# this only works if you are in Triton!!
from dmrgpy import parallelslurm
#ms = parallelslurm.pcall(getm,bs) # compute all the magnetic fields
ms = [getm(bi) for bi in bs] # the previous call is the parallelization of this

np.savetxt("M_VS_B.OUT",np.array([bs,ms]).T) # save in a file



