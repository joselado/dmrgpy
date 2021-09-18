# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')


import numpy as np
from dmrgpy import fermionchain

nf = 10 # number of different spinless fermionic orbitals
fc = fermionchain.Fermionic_Chain(nf) # create the object
C = fc.C # annihilation
Cdag = fc.Cdag # creation

H = 0 # initialize Hamiltonian 

# random Hamiltonian
for i in range(nf-1):
    H = H + Cdag[i]*C[i+1]*np.random.random() # random first neigh. hopping
    H = H + C[i]*C[i+1]*np.random.random() # random first neigh. pairing
    # random first neigh. interaction
    H = H + Cdag[i]*C[i]*Cdag[i+1]*C[i+1]*np.random.random() 

H = H + H.get_dagger() # make it Hermitian

fc.set_hamiltonian(H) # set the Hamiltonian
print("Energy with DMRG",fc.gs_energy(mode="DMRG")) # energy with DMRG
print("Energy with ED",fc.gs_energy(mode="ED")) # energy with exact diag.









