# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
import time
n = 16 # take n different sites
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Chain(spins) # create the spin chain
h = 0
for i in range(n-1):
    h = h + np.random.random()*sc.SS(i,i+1)
sc.set_hamiltonian(h)

nex = 3
sc.excited_gram_schmidt = True
es0,wfs0 = sc.get_excited_states(n=nex)
sc.excited_gram_schmidt = False
es1,wfs1 = sc.get_excited_states(n=nex)
# Hamiltonian times wavefunctions
hwfs0 = [h*wf for wf in wfs0]
hwfs1 = [h*wf for wf in wfs1]

de0 = [hwfs0[i].dot(hwfs0[i]) - wfs0[i].dot(hwfs0[i])**2 for i in range(nex)]
de1 = [hwfs1[i].dot(hwfs1[i]) - wfs1[i].dot(hwfs1[i])**2 for i in range(nex)]

de0 = np.sqrt(np.array(de0))
de1 = np.sqrt(np.array(de1))

print("Average fluctuation with GM",np.mean(de0))
print("Average fluctuation without GM",np.mean(de1))









