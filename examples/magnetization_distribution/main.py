# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain

n = 34 # number of sites in your chain
spins = [2 for i in range(n)] # create the sites
sc = spinchain.Spin_Chain(spins) # create the chain

# now define the Hamiltonian
h = 0
for i in range(n-1): h = h - sc.Sz[i]*sc.Sz[i+1] # add exchange
for i in range(n): h = h - .5*sc.Sx[i] # add transverse field
h = h + sc.Sz[n-1]*sc.Sz[0] # and apply periodic boundary conditions
sc.set_hamiltonian(h) # and initialize the Hamiltonian

# setup some parameters
sc.maxm = 60
sc.kpmmaxm = 40

# now define the operator for which you want the distribution
M = 0
for i in range(n): M += sc.Sz[i] # total magnetization 
x,y = sc.get_distribution(X=M,scale=n,delta=6e-1) # compute a distribution

# plot the result and save it in a file
import matplotlib.pyplot as plt
np.savetxt("DISTRIBUTION.OUT",np.array([x,y.real]).T)
plt.plot(x,y.real,marker="o") # correlator using DMRG
plt.show()


