# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
n = 5
fc = fermionchain.Fermionic_Chain(n) # create the chain
h = 0
#for i in range(n):
#  for j in range(i):
#    h = h + np.random.random()*fc.C[i]*fc.Cdag[j]
C = fc.A
Cdag = fc.Adag
for i in range(n):
  for j in range(n):
    for k in range(n):
      for l in range(n):
        Ni = C[i]*Cdag[j]
        Nj = C[k]*Cdag[l]
        h = h + np.random.random()*Ni*Nj
h = h + h.get_dagger()
#exit()

ms = [3,10,20,1000]
#ms = [10]
#fc.set_hamiltonian(h)
#e0 = fc.gs_energy(mode="ED") # energy with exact diagonalization
#print(e0)
import profile
#
#def f():
#  fc.set_hamiltonian(h)
#  e1 = fc.gs_energy(mode="DMRG") # energy with DMRG
#
#profile.run('f()','restats')
#import pstats
#from pstats import SortKey
#p = pstats.Stats('restats')
#p.sort_stats(SortKey.CUMULATIVE).print_stats(20)
#exit()
import time
fc.set_hamiltonian(h)
t0 = time.time()
e1 = fc.gs_energy(mode="DMRG") # energy with DMRG
t1 = time.time()
print("Time in GS = ",t1-t0)
#fc.mpomaxm = 500
fc.get_dynamical_correlator(name=(fc.C[0],fc.C[1]))
t2 = time.time()
print("Time in correlator = ",t2-t1)

