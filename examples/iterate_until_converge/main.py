# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
import fermionchain
n = 10
sc = fermionchain.Fermionic_Chain(n) # create the chain
def ft(i,j):
#    if i==j: return 1.0
    if abs(j-i)==1: return 1.0 #+ np.random.random()
    if i==j: return np.random.random()
#    if i==j: return 1.1
    return 0.0
sc.set_hoppings(ft) # hoppings
import time
#e0 = sc.gs_energy() # compute ground state energy with DMRG
e1 = sc.gs_energy_free() # compute ground state energy for free electrons
print("Free",e1)
sc.nsweeps = 1
es = [] # empty list
wf0 = None # no initial wavefunction
for i in range(10):
  sc.nsweeps = 1 # do nothing
  e = sc.gs_energy(wf0=sc.wf0) 
  es.append(e)
  print(e)
plt.plot(range(len(es)),es,marker="o")
plt.show()


