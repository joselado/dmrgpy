# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import fermionchain
n = 20 # number of spinful fermionic sites
fc = fermionchain.Spinful_Fermionic_Chain(n) # create the chain
h = 0
U = 2.0 # hubbard
J = 1.0 # local exchange

# define a linear hubbard chain, with a local magnetic field in an edge

for i in range(n-1): # first neighbor hopping
    h = h + fc.Cdagup[i]*fc.Cup[i+1]
    h = h + fc.Cdagdn[i]*fc.Cdn[i+1]
h = h + h.get_dagger() # add Hermitian conjugate
for i in range(n): # local Hubbard interaction
    h = h + U*(fc.Nup[i]-0.5)*(fc.Ndn[i]-0.5)

h = h + J*fc.Sz[0] # add the magnetic field in the edge

fc.set_hamiltonian(h) # set the Hamiltonian to the object

Mz = [fc.vev(Szi) for Szi in fc.Sz] # compute magnetization in the sites
inds = range(n) # indexes of the sites
import matplotlib.pyplot as plt
plt.plot(inds,Mz,marker="o")
plt.show()
