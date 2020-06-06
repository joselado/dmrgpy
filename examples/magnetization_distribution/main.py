# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain

n = 6 # number of sites in your chain
spins = [2 for i in range(n)] # create the sites
sc = spinchain.Spin_Chain(spins) # create the chain

# now define the Hamiltonian
h = 0
for i in range(n-1): h = h - sc.Sz[i]*sc.Sz[i+1] # add exchange
for i in range(n): h = h - .5*sc.Sx[i] # add transverse field
h = h - sc.Sz[n-1]*sc.Sz[0] # and apply periodic boundary conditions
sc.set_hamiltonian(h) # and initialize the Hamiltonian

# setup some parameters
sc.maxm = 60
sc.kpmmaxm = 40

#sc.itensor_version = "julia"

# now define the operator for which you want the distribution
M = 0
for i in range(n): M += sc.Sz[i] # total magnetization 
x,y = sc.get_distribution(X=M,delta=1e-1) # compute a distribution
x1,y1 = sc.get_distribution(X=M,delta=1e-1,mode="ED") # compute a distribution

# plot the result and save it in a file
print("Integral DMRG=",np.trapz(y.real,dx=x[1]-x[0]))
print("Integral ED=",np.trapz(y1.real,dx=x1[1]-x1[0]))
import matplotlib.pyplot as plt
np.savetxt("DISTRIBUTION.OUT",np.array([x,y.real]).T)
plt.plot(x,y.real,marker="o",label="DMRG") # correlator using DMRG
plt.plot(x1,y1.real,marker="o",label="ED") # correlator using DMRG
plt.xlabel("magnetization")
plt.ylabel("distribution")
plt.legend()
plt.show()


