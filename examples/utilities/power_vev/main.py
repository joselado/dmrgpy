# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain

n = 15 # number of sites in your chain
spins = ["S=1/2" for i in range(n)] # create the sites
sc = spinchain.Spin_Chain(spins) # create the chain
# Note: this example used to set sc.itensor_version = "julia", but that
# value refers to the legacy subprocess-based Julia path (juliarun.py),
# which is no longer reachable through the public API (only "julia_live"
# is) and raises inside is_hermitian()/random_mps() since self._session is
# never created for it. Use the default backend instead.


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

#wf = sc.get_gs()
#sc.applyoperator(h,wf)
#exit()

import time
t0 = time.time()
d1 = sc.vev(Mz,npow=2).real # compute Mz**2 with a more efficient algorithm
t1 = time.time()
d0 = sc.vev(Mz*Mz).real # compute Mz**2 by brute force
t2 = time.time()
print("Mz2 with smarter algorithm",d1,"and time spent",t1-t0)
print("Mz2 by brute force",d0,"and time spent",t2-t1)

print("\nPerformace improvement\n",(t2-t1)/(t1-t0))

# sweep the chain size and see how the performance improvement scales
def performance_improvement(ni):
    spins_i = ["S=1/2" for i in range(ni)]
    sci = spinchain.Spin_Chain(spins_i)
    hi = 0
    for i in range(ni-1):
        hi = hi - sci.Sz[i]*sci.Sz[i+1]
    for i in range(ni): hi = hi + 0.6*sci.Sx[i]
    Mzi = sum(sci.Sz)
    sci.set_hamiltonian(hi)
    sci.get_gs()
    t0 = time.time()
    sci.vev(Mzi,npow=2).real
    t1 = time.time()
    sci.vev(Mzi*Mzi).real
    t2 = time.time()
    return (t2-t1)/(t1-t0)

ns = range(4,16,3)
speedups = [performance_improvement(ni) for ni in ns]

import matplotlib.pyplot as plt

plt.plot(ns,speedups,marker="o")
plt.xlabel("Chain size n")
plt.ylabel("Performance improvement (brute force / smart)")
plt.show()











