# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
import time
n = 6 # take n different sites
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Chain(spins) # create the spin chain
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1]
    h = h + sc.Sy[i]*sc.Sy[i+1]
    h = h + sc.Sz[i]*sc.Sz[i+1]

# energies with ED
sc.set_hamiltonian(h)
es0 = sc.get_excited(n=6,mode="ED") # compute excited states with ED

for n in [4,6,8,10,20,30]: 
    sc.set_hamiltonian(h)
    sc.nsweeps = n
    scale = 10 # this is the lagrange multiplier, default is 10
    es = sc.get_excited(n=6,scale=scale) # compute excited states with DMRG
    error = np.mean(np.abs(np.array(es0) - np.array(es))) # error
    print(n,error)
