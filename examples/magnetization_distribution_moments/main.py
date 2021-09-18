# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain

n = 20 # number of sites in your chain
spins = ["S=1/2" for i in range(n)] # create the sites
sc = spinchain.Spin_Chain(spins) # create the chain
#sc.itensor_version = "julia"

# now define the Hamiltonian for a transverse Ising model
h = 0 # initialization
for i in range(n-1): 
    h = h - sc.Sz[i]*sc.Sz[i+1] # add exchange
for i in range(n): 
    h = h - .3*sc.Sx[i] # add transverse field
sc.set_hamiltonian(h) # initialize the Hamiltonian
# setup some parameters
sc.maxm = 20 # bond dimension for ground state
sc.kpmmaxm = sc.maxm # bond dimension for KPM algorithm


# now define the operator for which you want the distribution
M = (1/n)*sum(sc.Sz) # normalized magnetization

# you can just compute the distribution with
#x,y = sc.get_distribution(X=M,delta=2e-2) # compute a distribution


# compute the moments of the distribution
mus,shift,scale = sc.get_distribution_moments(X=M,delta=2e-2) 
# delta controls the number of polynomials, make it bigger
# to compute less polynomials

# once we have the moments, lets reconstruct the function
# the available kernels are jackson, lorentz and plain
from dmrgpy.algebra.kpm import reconstruct_chebyshev
x,y = reconstruct_chebyshev(mus,shift=shift,
        scale=scale,x=np.linspace(-.7,.7,2000),kernel="jackson")

import matplotlib.pyplot as plt
np.savetxt("DISTRIBUTION.OUT",np.array([x,y.real]).T)
plt.plot(x,y.real,marker="o") # distribution using DMRG
plt.xlabel("magnetization")
plt.ylabel("distribution")
plt.show()











