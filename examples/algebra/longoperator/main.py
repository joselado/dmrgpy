# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import fermionchain
ns = 20 # number of spinful fermionic sites
fc = fermionchain.Spinful_Fermionic_Chain(ns) # create the chain
h = 0 # initialize Hamiltonian

mu = -1.0 # chemical potential

# first add the hoppings
for i in range(ns-1): # loop over sites
    h = h + fc.Cdagup[i]*fc.Cup[i+1]
    h = h + fc.Cdagdn[i]*fc.Cdn[i+1]
h = h + h.get_dagger() # complex conjugate


# now add the chemical potential
for i in range(ns): # loop over sites
    h = h + mu*fc.Cdagup[i]*fc.Cup[i]
    h = h + mu*fc.Cdagdn[i]*fc.Cdn[i]

# setup the Hamiltonian
fc.set_hamiltonian(h) # initialize Hamiltonian
wf = fc.get_gs() # get ground state
i = 0
inds = range(ns)

# function to enable computing correlators in big systems
from dmrgpy.multioperatortk.longoperator import toMPO

vevs0,vevs1 = [],[]
for j in inds: # loop over correlators
    op0 = fc.Cdagup[i]*fc.Cup[j] # operator
    op1 = toMPO(fc,op0) # use the prodcut trick (enables very big systems)
    vevs0.append(wf.dot(op0*wf)) # full mode
    vevs1.append(wf.dot(op1*wf)) # prodcut mode


import matplotlib.pyplot as plt
plt.plot(inds,vevs1,marker="o",label="product")
plt.plot(inds,vevs0,label="full")
plt.legend()
plt.xlabel("Site j")
plt.ylabel("$\\langle C^\\dagger_i C_j \\rangle$")
plt.show()

