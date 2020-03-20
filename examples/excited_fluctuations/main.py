# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
import time
n = 16 # take n different sites
spins = [2 for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Chain(spins) # create the spin chain
h = 0
for i in range(n-1):
    h = h + np.random.random()*sc.SS(i,i+1)
sc.set_hamiltonian(h)

sc.excited_gram_schmidt = True
de0 = sc.get_excited(n=20,fluctuations=True)[1]
sc.excited_gram_schmidt = False
de1 = sc.get_excited(n=20,fluctuations=True)[1]

print("Average fluctuation with GM",np.mean(de0))
print("Average fluctuation without GM",np.mean(de1))
