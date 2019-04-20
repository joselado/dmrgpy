# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
n = 10
fc = fermionchain.Fermionic_Hamiltonian(n,spinful=False) # create the chain
m = np.matrix(np.random.random((n,n)) + 1j*np.random.random((n,n)))
m = m + m.H
def ft(i,j):
    return m[i,j]
    if abs(j-i)==1: return 1.0 
    return 0.0
fc.set_hoppings(ft) # hoppings
pairs = [(0,i) for i in range(n)]
out1 = fc.get_correlator(pairs=pairs,mode="DMRG") # using the default function
pvev = [[("Cdag",p[0]),("C",p[1]) ] for p in pairs] # define operators
out2 = np.array([fc.vev(p) for p in pvev]) # compute using the VEV feature

import matplotlib.pyplot as plt

plt.plot(range(len(pairs)),out1.real,c="red")
plt.scatter(range(len(pairs)),out2.real,c="blue")
plt.show()


