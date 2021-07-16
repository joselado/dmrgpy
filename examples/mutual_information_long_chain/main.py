# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
n = 20
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain

# This example shows how to compute the mutual information

# Heisenberg model
sc = spinchain.Spin_Chain(spins) # create the spin chain
h = 0
for i in range(n-1):
    h = h +sc.Sx[i]*sc.Sx[i+1]
    h = h +sc.Sy[i]*sc.Sy[i+1]
    h = h +sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
wf = sc.get_gs() # compute ground state

inds = range(n-1)
ss = [wf.get_mutual_information(i,i+1) for i in inds]

import matplotlib.pyplot as plt

plt.plot(inds,ss,marker="o")
plt.xlabel("Bond")
plt.ylabel("Mutual information")
plt.show()

