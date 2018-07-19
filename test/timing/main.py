from __future__ import print_function
import sys
import os
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg
import spinchain
import time

ns = range(10,14,2)
ts = []
for n in ns:
#  spins = [np.random.randint(2,7) for i in range(n)] # spin 1/2 heisenberg chain
  spins = [2 for i in range(n)] # spin 1/2 heisenberg chain
  sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain
  t0 = time.time()
#  e0 = sc.gs_energy(mode="DMRG") # compute the ground state energy
#  e0 = sc.get_excited(n=4,mode="DMRG") # compute the ground state energy
  sc.get_dynamical_correlator(mode="DMRG",i=1,j=3,name="XX",delta=0.02)
  t1 = time.time()
  ts.append(t1-t0) # store
  print("Time in ",n,t1-t0)
#  e1 = sc.gs_energy(mode="ED") # compute the ground state energy
#  print("Spin chain",spins)
#  print("Energy with ED",e1)
#  print("Energy with DMRG",e0)
#  print("\n")
#  de = np.abs(e1-e0)
#  if de>0.0001: raise

#print("Test passed")
ts = np.array(ts)
ns = np.array(ns)



out = np.polyfit(np.log(ns),np.log(ts),deg=1)

print("Power",out[0])

import matplotlib.pyplot as plt

plt.plot(np.log(ns),np.log(ts),marker="o")
plt.xlabel("log(N)")
plt.ylabel("log(T)")
plt.title("Power ="+str(out[0]))

plt.show()
