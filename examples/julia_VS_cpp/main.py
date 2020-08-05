# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain

def get(version,n=30):
  spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
  sc = spinchain.Spin_Chain(spins) # create the spin chain
  if version=="julia": sc.setup_julia() # setup the Julia mode
  h = 0
  for i in range(n-1):
      h = h +sc.Sx[i]*sc.Sx[i+1]
      h = h +sc.Sy[i]*sc.Sy[i+1]
      h = h +sc.Sz[i]*sc.Sz[i+1]
  sc.set_hamiltonian(h)
  return sc.gs_energy()


e0 = get("julia",n=10) # warmup round with Julia to precompile

import time

def time_ratio(n):
  t0 = time.time()
  e0 = get(2,n=n)
  t1 = time.time()
  e1 = get("julia",n=n)
  t2 = time.time()
  print("Energy with C++",e0)
  print("Energy with Julia",e1)
  print("Time with C++",t1-t0)
  print("Time with Julia",t2-t1)
  return (t1-t0)/(t2-t1)


for n in range(10,30,2):
    print("Length",n,"C++/Julia",time_ratio(n))


