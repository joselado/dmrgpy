# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import fermionchain
n = 20 # number of fermionic sites
fc = fermionchain.Fermionic_Chain(n) # create the chain
h = 0
U = -0.95
for i in range(n-1): # hopping
    h = h + fc.Cdag[i]*fc.C[i+1]
    h = h + U*(fc.N[i]-.5)*(fc.N[i+1]-.5)
h = h + h.get_dagger()
##############################
# Setup the Many Body Hamiltonian
fc.maxm = 20
fc.nsweeps = 10
fc.set_hamiltonian(h) # set the hoppings
wf = fc.get_gs(mode="DMRG")

#from dmrgpy.mps import random_mps
#wf = random_mps(fc)

import time

t0 = time.time()
m1 = wf.get_correlation_matrix(dmmode="fast")
t1 = time.time()


m2 = wf.get_correlation_matrix(dmmode="full")
t2 = time.time()

print("Python mode",t1-t0)
print("C++ mode",t2-t1)

print("Difference",np.mean(np.abs(m1-m2)))

