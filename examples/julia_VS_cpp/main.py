# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
n = 20
spins = [2 for i in range(n)] # spin 1/2 heisenberg chain
def get():
  sc = spinchain.Spin_Chain(spins) # create the spin chain
  h = 0
  for i in range(n-1):
      h = h +sc.Sx[i]*sc.Sx[i+1]
      h = h +sc.Sy[i]*sc.Sy[i+1]
      h = h +sc.Sz[i]*sc.Sz[i+1]
  sc.set_hamiltonian(h)
  return sc


import time
t0 = time.time()
sc = get()
e0 = sc.gs_energy() # compute the ground state energy
t1 = time.time()
sc = get()
sc.itensor_version = "julia" # setup this version
e1 = sc.gs_energy() # compute the ground state energy
t2 = time.time()
print("Energy with C++",e0)
print("Energy with Julia",e1)
print("Time with C++",t1-t0)
print("Time with Julia",t2-t1)
