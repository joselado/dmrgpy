# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import fermionchain
n = 4
# This example shows how compute the entropy between sites (0,1) and

# function computing the entropy
def get(c01):
    """Compute entropy"""
    fc = fermionchain.Fermionic_Chain(n)
    h = 0
    for i in range(n-1):
        c = 1.0 # default coupling
        h = h -fc.C[i]*fc.Cdag[i]
    h = h + h.get_dagger() # add the Hermitian conjugate
    fc.set_hamiltonian(h)
    wf = fc.get_gs() # compute ground state
    s = wf.get_pair_entropy(0,None) # entropy of the sites (0,1) with the rest
    return s

print(get(0.)) ; exit()

cs = np.linspace(0.,1.,20) # coupling between left and right parts
ss = [get(c) for c in cs] # compute all the entropies

import matplotlib.pyplot as plt

plt.plot(cs,ss,marker="o")
plt.xlabel("Coupling between (1,2) and (3,4)")
plt.ylabel("Entanglement entropy of the (1,2) pair")
plt.show()










