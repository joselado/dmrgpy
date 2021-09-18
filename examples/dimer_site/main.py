# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import fermionchain
n = 2 # number of spinful fermionic sites
fc = fermionchain.Fermionic_Chain(n) # create the chain
h = fc.N[0] + fc.N[1]
fc.set_hamiltonian(h)
##############################
# Setup the Many Body Hamiltonian
print(fc.gs_energy(mode="DMRG"))
print(fc.gs_energy(mode="ED"))
print("N",fc.vev(fc.N[0],mode="DMRG"))
print("N",fc.vev(fc.N[0],mode="ED"))
j = 0 # second index of the dynamical correlator
delta = 0.1 # energy resolution (approximate)
name = (fc.C[0],fc.Cdag[0])
(x,y) = fc.get_dynamical_correlator(delta=delta,name=name)
(x1,y1) = fc.get_dynamical_correlator(delta=delta,name=name,mode="ED",
        submode="INV")
import matplotlib.pyplot as plt
# plot the result
plt.plot(x,y.real,marker="o",label="DMRG")
plt.plot(x1,y1.real,marker="o",label="ED")
plt.legend()
plt.show()











