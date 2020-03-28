# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import fermionchain
n = 1 # number of spinful fermionic sites
fc = fermionchain.Spinful_Fermionic_Chain(n) # create the chain
h = fc.Cup[0]*fc.Cdn[0]
h = h + h.get_dagger()
fc.set_hamiltonian(h)
print(fc.get_excited(n=2,mode="ED"))
##############################
# Setup the Many Body Hamiltonian
print(fc.gs_energy(mode="DMRG"))
print(fc.gs_energy(mode="ED"))
print("Sz",fc.vev(fc.Sz[0],mode="DMRG"))
print("Sz",fc.vev(fc.Sz[0],mode="ED"))
j = 0 # second index of the dynamical correlator
delta = 0.1 # energy resolution (approximate)
name = (fc.Cup[0],fc.Cdagup[0])
(x,y) = fc.get_dynamical_correlator(delta=delta,name=name)
(x1,y1) = fc.get_dynamical_correlator(delta=delta,name=name,mode="ED",
        submode="INV")
import matplotlib.pyplot as plt
# plot the result
plt.plot(x,y.real,marker="o",label="DMRG")
plt.plot(x1,y1.real,marker="o",label="ED")
plt.legend()
plt.show()


