import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain

nf = 8 # number of different spinless fermionic orbitals
fc = fermionchain.Fermionic_Hamiltonian(nf) # create the object
C = fc.C # annihilation
Cdag = fc.Cdag # creation
N = fc.N # density
#N = [Cdag[i]*C[i] for i in range(nf)] # density

H = 0 # initialize Hamiltonian 

# random Hamiltonian
for i in range(nf-1):
    H = H + Cdag[i]*C[i+1]*np.random.random() # random first neigh. hopping
    H = H - N[i]*N[i+1]*np.random.random() # random first neigh. interaction

H = H + H.get_dagger() # make it Hermitian

fc.set_hamiltonian(H)
print("Energy with DMRG",fc.gs_energy(mode="DMRG"))
print("Energy with ED",fc.gs_energy(mode="ED"))
