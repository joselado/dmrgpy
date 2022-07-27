# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
n = 8
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Chain(spins) # create the spin chain

# fidelity of the transverse field Ising model
def get(l):
    # original Hamiltonian
    h0 = 0
    for i in range(n-1):
        h0 = h0 + 4*sc.Sz[i]*sc.Sz[i+1]
    h0 = h0 + 4*sc.Sz[0]*sc.Sz[n-1]
    # perturbation
    h1 = 0
    print(l)
    for i in range(n): h1 = h1 + 2*sc.Sx[i]
    from dmrgpy import fidelity
    return fidelity.get_fidelity(sc,h0,h1,l,n=3,mode="ED")

#get(0.5) ; exit()

ls = np.linspace(0.,2.,10)
fs = [get(l) for l in ls]

import matplotlib.pyplot as plt
plt.xlabel("Bx")
plt.ylabel("Fidelity susceptibility")
plt.plot(ls,fs)
plt.show()








