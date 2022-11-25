# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import fermionchain
def get(gamma):
    n = 8
    fc = fermionchain.Fermionic_Chain(n) # create the fermion chain
    mh = np.zeros((n,n),dtype=np.complex) # TB matrix
#    gamma = 0.2
    V = -0.3
    for i in range(n-1):
        mh[i,i+1] = 1.0 + gamma
        mh[i+1,i] = 1.0 - gamma
    h = 0 # initialize Hamiltonian
    for i in range(n-1): 
        h = h + V*(fc.N[i]-0.5)*(fc.N[i+1]-0.5)
    for i in range(n):
        for j in range(n):
            h = h + mh[i,j]*fc.Cdag[i]*fc.C[j]
    fc.set_hamiltonian(h) 
    es = fc.get_excited(mode="ED",n=10)
    return es-es[0].real

import matplotlib.pyplot as plt

gs = np.linspace(0.,1.5,20)
for g in gs:
    es = get(g)
    plt.subplot(1,2,1)
    plt.scatter(es*0+g,es.real,c="black")
    plt.subplot(1,2,2)
    plt.scatter(es*0+g,es.imag,c="black")
plt.show()








