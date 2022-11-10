# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain

def get(delta):
    n = 4
    spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
    sc = spinchain.Spin_Chain(spins) # create the spin chain
    h = 0
    for i in range(n-1):
        h = h +sc.Sx[i]*sc.Sx[i+1]
        h = h +sc.Sy[i]*sc.Sy[i+1]
        h = h +sc.Sz[i]*sc.Sz[i+1]
    
    for i in range(n):
        h = h + (-1)**i*sc.Sz[i]*0.6j # add some imaginary part
    sc.set_hamiltonian(h) # set Hamiltonian
    return sc.get_gs_degeneracy(mode="ED",delta=delta)

deltas = np.linspace(1e-5,0.4,30)
degs = [get(delta) for delta in deltas]

import matplotlib.pyplot as plt

plt.plot(deltas,degs)
plt.xlabel("Energy smearing")
plt.ylabel("Degeneracy")
plt.show()


