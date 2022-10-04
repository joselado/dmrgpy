# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain
n = 6
# create a random spin chain
#spins = [np.random.randint(2,5) for i in range(n)] # spin 1/2 heisenberg chain
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain


# create first neighbor exchange
sc = spinchain.Spin_Chain(spins) # create the spin chain
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1]
    h = h + sc.Sy[i]*sc.Sy[i+1]
    h = h + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)

#print(sc.gs_energy())

sc.kpmmaxm = 20 # bond dimension for KPM
sc.maxm = 20 # bond dimension

name = (sc.Sz[0],sc.Sz[0]) # correlator to compute
es = np.linspace(-0.5,4,200)
delta = 5e-2
(x,y) = sc.get_dynamical_correlator(name=name,es=es,delta=delta)

import matplotlib.pyplot as plt
plt.plot(x,y.real,c="blue",label="DMRG")
plt.show()












