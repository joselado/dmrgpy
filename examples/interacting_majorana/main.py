# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import fermionchain
n = 20 # number of spinful fermionic sites
fc = fermionchain.Fermionic_Chain(n) # create the chain
h = 0
V = -0.9
for i in range(n-1): # hopping
    h = h + fc.Cdag[i]*fc.C[i+1]
#    h = h + 0.5*(-1)**i*fc.Cdag[i]*fc.C[i+1]
#    h = h + 0.3*fc.C[i]*fc.C[i+1]
    h = h + V*(fc.N[i]-0.5)*(fc.N[i+1]-0.5)
h = h + h.get_dagger()
##############################
# Setup the Many Body Hamiltonian
fc.maxm = 20
fc.kpmmaxm = fc.maxm
fc.nsweeps = 10
fc.set_hamiltonian(h) # set the hoppings
fc.get_gs()

ii = 0#n//2
(es,ds) = fc.get_dynamical_correlator(es=np.linspace(0.,3.,600),name=(fc.Cdag[ii],fc.C[ii]),delta=1e-1)

np.savetxt("DOS.OUT",np.array([es,ds.real]).T)









