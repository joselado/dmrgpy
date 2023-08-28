# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import fermionchain
ns = 10 # number of spinful fermionic sites
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

i = 0
inds = range(ns)
vevup = [fc.vev(fc.Cdagup[i]*fc.Cup[j]) for j in inds]
vevdn = [fc.vev(fc.Cdagdn[i]*fc.Cdn[j]) for j in inds]


# these lines would be for charge-charge correlator
#vevup = [fc.vev(fc.Nup[i]*fc.Nup[j]) - fc.vev(fc.Nup[i])*fc.vev(fc.Nup[j]) for j in inds]
#vevdn = [fc.vev(fc.Ndn[i]*fc.Ndn[j]) - fc.vev(fc.Ndn[i])*fc.vev(fc.Ndn[j]) for j in inds]




vevs = np.array(vevup) + np.array(vevdn) # up + down



import matplotlib.pyplot as plt
plt.plot(inds,vevs,marker="o")
plt.xlabel("Site j")
plt.ylabel("$\\langle C^\\dagger_i C_j \\rangle$")
plt.show()

