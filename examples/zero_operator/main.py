# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import parafermionchain
n = 4
pc = parafermionchain.Parafermionic_Chain(n) # create the chain

# define an operator that should be zero due to the anticommutation relations
for ii in range(ntries):
    i = np.random.randint(n) # get one index
    j = np.random.randint(n) # get another index

h = h +h.get_dagger() # Make the Hamiltonian Hermitian
pc.set_hamiltonian(h) # set the Hamiltonian
print(pc.gs_energy(mode="DMRG"))
print(pc.gs_energy(mode="ED"))
