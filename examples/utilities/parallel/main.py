# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain

def get(u):
    """Ground state energy of a Heisenberg chain with a staggered
    field of amplitude u, so that each parallel call solves a
    physically distinct Hamiltonian"""
    n = 16
    spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
    sc = spinchain.Spin_Chain(spins) # create the spin chain
    h = 0
    for i in range(n-1):
        h = h +sc.Sx[i]*sc.Sx[i+1]
        h = h +sc.Sy[i]*sc.Sy[i+1]
        h = h +sc.Sz[i]*sc.Sz[i+1]
    for i in range(n):
        h = h + u*((-1)**i)*sc.Sz[i] # staggered field, distinct per worker

    sc.set_hamiltonian(h)
    from dmrgpy import algebra ; algebra.maxsize = 70000
    #print(sc.get_excited(n=1,mode="ED")) ; exit()
    e = sc.gs_energy() # compute the ground state energy
    print("Energy",e)
    return e

from dmrgpy import parallel
parallel.cores = 4
us = np.linspace(0.,1.,4) # staggered field amplitude for each parallel call
es = parallel.pcall(get,us)

import matplotlib.pyplot as plt

plt.plot(us,es,marker="o")
plt.xlabel("Staggered field u")
plt.ylabel("Ground state energy")
plt.show()









