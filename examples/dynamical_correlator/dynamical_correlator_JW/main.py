# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain,fermionchain
n = 6

site = 0

def get_spin():
    # create a random spin chain
    spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
    # create first neighbor exchange
    sc = spinchain.Spin_Chain(spins) # create the spin chain
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1]
        h = h + sc.Sy[i]*sc.Sy[i+1]
    sc.set_hamiltonian(h)
    name = (sc.Sx[site],sc.Sx[site]) # correlator to compute
    es = np.linspace(-0.5,3,2000)
    delta = 5e-2
    (x,y) = sc.get_dynamical_correlator(name=name,es=es,delta=delta)
    return (x,y)


def get_fermion():
    sc = fermionchain.Fermionic_Chain(n) # create the spin chain
    h = 0
    for i in range(n-1):
        h = h + sc.Cdag[i]*sc.C[i+1]/2.
    h = h + h.get_dagger()
    sc.set_hamiltonian(h)
    name = (sc.Cdag[site],sc.C[site]) # correlator to compute
    es = np.linspace(-0.5,3,2000)
    delta = 5e-2
    (x,y) = sc.get_dynamical_correlator(name=name,es=es,delta=delta)
    return (x,y/2)

(xs,ys) = get_spin()
(xf,yf) = get_fermion()


import matplotlib.pyplot as plt
plt.plot(xf,yf.real,c="blue",label="Fermion")
plt.scatter(xs,ys.real,c="red",label="Spin")
plt.legend()
plt.show()












