# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import fermionchain
n = 4 # number of spinful fermionic sites
fc = fermionchain.Spinful_Fermionic_Chain(n) # create the chain
#fc.Nup = fc.N
#fc.Ndn = fc.N
#fc.Cup = fc.C
#fc.Cdagup = fc.Cdag
#fc.Cdn = fc.C
#fc.Cdagdn = fc.Cdag
# initialize Hamiltonian #
h = 0
for i in range(n-1): # hopping
    h = h + fc.Cdagup[i]*fc.Cdagup[i+1]
    h = h + fc.Cdn[i]*fc.Cdagdn[i+1]
for i in range(n): # Hubbard
    h = h + (fc.Nup[i]-.5)*(fc.Ndn[i]-.5)
h = h + h.get_dagger()
##############################
# Setup the Many Body Hamiltonian
fc.maxm = 40
fc.nsweeps = 40
fc.set_hamiltonian(h) # set the hoppings
print(fc.gs_energy(mode="DMRG"))
print(fc.gs_energy(mode="ED"))
#exit()
# Compute the dynamical correlator defined by
# <0|c_i^dagger \delta(H-E_0-\omega) c_j |0>
i = 0 # first index of the dynamical correlator
j = 0 # second index of the dynamical correlator
delta = 0.1 # energy resolution (approximate)
fc.kpmmaxm = 20 # maximum bond dimension in KPM
fc.gs_energy()
#exit()
# The result will be written in a file called DYNAMICAL_CORRELATOR.OUT
# compute the dynamical correlator using KPM DMRG
name = (fc.N[0],fc.N[0])
(x,y) = fc.get_dynamical_correlator(delta=delta,name=name)
import matplotlib.pyplot as plt
# plot the result
plt.plot(x,y.real,marker="o")
plt.show()


