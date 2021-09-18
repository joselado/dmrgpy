# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain
n = 4
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain

# This example shows how compute the entropy between sites (0,1) and
# (2,3) of a spin chain with 4 sites

# function computing the entropy
def get(c01):
    """Compute entropy"""
    # Heisenberg model
    sc = spinchain.Spin_Chain(spins) # create the spin chain
    h = 0
    for i in range(n-1):
        c = 1.0 # default coupling
        if i==1: c = c01 # if it is the middle, change the coupling
        h = h +c*sc.Sx[i]*sc.Sx[i+1]
        h = h +c*sc.Sy[i]*sc.Sy[i+1]
        h = h +c*sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    wf = sc.get_gs() # compute ground state
    s = wf.get_pair_entropy(0,1) # entropy of the sites (0,1) with the rest
#   Alternative method
#    s1 = wf.get_bond_entropy(1,2) # compute entropy at this bond
#    print(s,s1,s1/s)
    return s

cs = np.linspace(0.,1.,20) # coupling between left and right parts
ss = [get(c) for c in cs] # compute all the entropies

import matplotlib.pyplot as plt

plt.plot(cs,ss,marker="o")
plt.xlabel("Coupling between (1,2) and (3,4)")
plt.ylabel("Entanglement entropy of the (1,2) pair")
plt.show()










