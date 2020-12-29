# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
import time
n = 8 # take n different sites
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Chain(spins) # create the spin chain
h = 0
for i in range(n-1):
    h = h + sc.SS(i,i+1)

def get_ex(mode="DMRG",purify=True,nsweeps=10):
    sc.maxm = 10
    sc.nsweeps = nsweeps
    sc.set_hamiltonian(h)
    if mode=="DMRG":
        return sc.get_excited(n=10,mode="DMRG",purify=purify) # compute ex. st.
    if mode=="ED":
        return sc.get_excited(n=10,mode="ED") # compute ex. st.


esed = get_ex(mode="ED") # get the "true" excited states with ED
diffpure = [] # empty list
diff = [] # empty list
ns = range(2,7,2) # different number of sweeps

for ii in ns: # loop over sweeps
  es = get_ex(mode="DMRG",purify=False,nsweeps=ii) # no purification
  diff.append(np.mean(np.abs(es-esed))) # error with respect to ED
  es = get_ex(mode="DMRG",purify=True,nsweeps=ii) # no purification
  diffpure.append(np.mean(np.abs(es-esed))) # error with respect to ED

import matplotlib.pyplot as plt

plt.plot(ns,diff,label="No purification",marker="o")
plt.plot(ns,diffpure,label="Purification",marker="o")
plt.xlabel("Sweeps")
plt.ylabel("log(error)")
plt.legend()
plt.show()






