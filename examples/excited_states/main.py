# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import spinchain
import time
n = 12 # take n different sites
spins = [2 for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain
# n is the number of excited states
t0 = time.time()
es1 = sc.get_excited(n=6,mode="DMRG") # compute excited states with DMRG
t1 = time.time()
es2 = sc.get_excited(n=6,mode="ED") # compute excited states with ED
t2 = time.time()
print("Time spent in DMRG",t1-t0)
print("Excited states with DMRG")
print(es1)
print("\n")
print("Time spent in ED",t2-t1)
print("Excited states with ED")
print(es2)
print("\n")
print("Relative error",np.sum(np.abs(es1-es2)))
sc.clean() # remove temporal files


