# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain


def get_dc(bc="open",n=3):
    spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
    sc = spinchain.Spin_Chain(spins) # create the spin chain
    h = 0
    if bc == "open": NN = n-1
    elif bc == "closed": NN = n
    print(NN)
    for i in range(NN):
        h = h + sc.Sx[i]*sc.Sx[(i+1)%n]
        h = h + sc.Sy[i]*sc.Sy[(i+1)%n]
        h = h + sc.Sz[i]*sc.Sz[(i+1)%n]
    Bz = 0.
    h = h + Bz*sum(sc.Sz)
    sc.set_hamiltonian(h)
    
    es = np.linspace(-0.5,4,200)
    delta = 1e-1
    (x,y) = sc.get_full_SS_correlator(es=es,delta=delta)
    return (x,y)

import matplotlib.pyplot as plt
bcs = ["open","closed"]
ns = [3,5]
ii = 0 
fig = plt.figure(figsize=(6,4))
for n in ns: # loop over system sizes
    plt.subplot(1,len(ns),ii+1) ; ii+=1
    plt.title("Size "+str(n))
    for bc in bcs: # loop over boundary conditions
        (x,y) = get_dc(bc=bc,n=n)
        plt.plot(x,y.real,label=bc,linewidth=3)
    plt.legend()
    plt.xlabel("Energy") ; plt.ylabel("Spectral fucntion")

plt.tight_layout()
plt.show()












