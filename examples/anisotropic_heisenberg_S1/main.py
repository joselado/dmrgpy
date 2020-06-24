# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
n = 40
spins = ["S=1/2"]+["S=1" for i in range(n)]+["S=1/2"] # spin 1/2 heisenberg chain
spins = ["S=1" for i in range(n)] # spin 1 heisenberg chain

def get(D):
    sc = spinchain.Spin_Chain(spins) # create the spin chain
    h = 0
    for i in range(len(spins)-1):
        h = h +sc.Sx[i]*sc.Sx[i+1]
        h = h +sc.Sy[i]*sc.Sy[i+1]
        h = h +sc.Sz[i]*sc.Sz[i+1]
    for i in range(len(spins)):  h = h - D*sc.Sz[i]*sc.Sz[i]
    sc.set_hamiltonian(h)
    sc.nsweeps = 5
    s = sc.get_gs().get_entropy(len(spins)//2)
    return s



import matplotlib.pyplot as plt


for d in np.linspace(-2.0,2.0,40):
    print(d)
    s = get(d)
    plt.scatter(s*0.+d,s,c="black")

plt.show()
