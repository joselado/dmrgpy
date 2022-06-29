# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain
n = 2
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain


sc = spinchain.Spin_Chain(spins) # create the spin chain

# create the Hamiltonian
Sx,Sy,Sz = sc.Sx,sc.Sy,sc.Sz
# add Heisenberg coupling between two sites
def SS(i,j): return  Sx[i]*Sx[j] + Sy[i]*Sy[j] + Sz[i]*Sz[j]

J = 1.0 # coupling between sites 0 and 1
h = SS(0,1) # add exchange 
sc.set_hamiltonian(h)

name = (Sz[0],Sz[0]) # correlator to compute
es = np.linspace(-0.5,4,200) # energy grid
delta = 5e-2 # smearing
(es,dc) = sc.get_dynamical_correlator(name=name,es=es,delta=delta)

import matplotlib.pyplot as plt
plt.plot(es,dc.real,c="blue",label="DMRG")
plt.show()












