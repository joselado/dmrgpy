# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import fermionchain
n = 4 # number of fermionic sites
fc = fermionchain.Fermionic_Chain(n) # create the chain
h = 0
U = 1.
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

fc.set_hamiltonian(h) # set the hoppings
wfed = fc.get_gs(mode="ED")

#from dmrgpy.mps import random_mps
#wf = random_mps(fc)

import time

# the correlation tensor is computed as
# M = <Cdag_i C_j Cdag_k C_l>

# first the full explicit MPS Python version
t0 = time.time()
m1 = wf.get_four_correlation_tensor(ctmode="explicit")
t1 = time.time()


# now the accelerated MPS algorithm
m2 = wf.get_four_correlation_tensor(ctmode="full")
t2 = time.time()

# ED as a reference

med = wfed.get_four_correlation_tensor()


# print the times
print("Python mode",t1-t0)
print("C++ mode",t2-t1)

def diff(m1,m2): return np.round(np.mean(np.abs(m1-m2)),5)

print("Difference full MPS and explicit MPS",diff(m1,m2))
print("Difference full MPS and ED",diff(m2,med))

