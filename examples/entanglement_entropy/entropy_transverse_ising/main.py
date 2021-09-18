# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain

def gets(b):
    n = 30
    spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
    sc = spinchain.Spin_Chain(spins) # create the spin chain
    h = 0
    for i in range(n-1): h = h - sc.Sz[i]*sc.Sz[i+1]
    for i in range(n): h = h - b*sc.Sx[i]
    
    sc.set_hamiltonian(h)
    wf = sc.get_gs() # compute ground state
    return wf.get_bond_entropy(n//2)

import matplotlib.pyplot as plt
bs = np.linspace(0.,1.0,30)
ss = [gets(b) for b in bs]
plt.plot(bs,ss,marker="o")
plt.show()










