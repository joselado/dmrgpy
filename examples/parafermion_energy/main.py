# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import parafermionchain
n = 4
pc = parafermionchain.Parafermionic_Chain(n) # create the chain

h = 0
for i in range(n):
  for j in range(n):
      h = h + pc.N[i]*pc.N[j]*np.random.random()
      h = h + pc.Sig[i]*pc.Sig[j]*np.random.random()
      h = h + pc.Sig[i]*pc.Sigd[j]*np.random.random()
      h = h + pc.Tau[i]*pc.Tau[j]*np.random.random()
      h = h + pc.Tau[i]*pc.Taud[j]*np.random.random()

h = h +h.get_dagger() # Make the Hamiltonian Hermitian
pc.set_hamiltonian(h) # set the Hamiltonian
print(pc.gs_energy(mode="DMRG"))
print(pc.gs_energy(mode="ED"))
