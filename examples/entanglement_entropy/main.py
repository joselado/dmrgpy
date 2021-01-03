# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
n = 10
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Chain(spins) # create the spin chain
h = 0
for i in range(n-1):
    h = h +sc.Sx[i]*sc.Sx[i+1]
    h = h +sc.Sy[i]*sc.Sy[i+1]
    h = h +sc.Sz[i]*sc.Sz[i+1]

# Add a magnetic field in the first site to polarize it
#h = h + 5.*sc.Sz[0]

sc.set_hamiltonian(h)
wf = sc.get_gs() # compute ground state
# Compute the entanglement entropy of every site
#ss = [wf.get_site_entropy(i) for i in range(n)]
ss = [wf.get_bond_entropy(0,i) for i in range(1,n)]

import matplotlib.pyplot as plt

plt.plot(range(len(ss)),ss,marker="o")
plt.show()

