# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain

n = 40 # number of sites in your chain
spins = ["S=1/2" for i in range(n)] # create the sites
sc = spinchain.Spin_Chain(spins) # create the chain

# now define the Hamiltonian
def geth(b):
  h = 0
  for i in range(n-1): 
      h = h - sc.Sz[i]*sc.Sz[i+1] # add exchange
  for i in range(n): h = h + b*sc.Sx[i] # add transverse field
  return h

# define the total spin in the z-direction
Mz = 0
for i in range(n): Mz = Mz + sc.Sz[i]

h = geth(0.6) # get the Hamiltonian
sc.set_hamiltonian(h) # and initialize the Hamiltonian
sc.get_gs() # compute the ground state

import time
t0 = time.time()
d1 = sc.vev(Mz,npow=2).real # compute Mz**2 with a more efficient algorithm
t1 = time.time()
d0 = sc.vev(Mz*Mz).real # compute Mz**2 by brute force
t2 = time.time()
print("Mz2 with smarter algorithm",d1,"and time spent",t1-t0)
print("Mz2 by brute force",d0,"and time spent",t2-t1)

print("\nPerformace improvement\n",(t2-t1)/(t1-t0))


