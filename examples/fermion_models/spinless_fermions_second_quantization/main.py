# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')


import numpy as np
from dmrgpy import fermionchain

def get_energies(nf,seed):
    """Ground-state energy (DMRG and ED) of a random spinless-fermion
    Hamiltonian on nf orbitals, for a fixed random seed (reproducible)."""
    np.random.seed(seed)
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
    edmrg = fc.gs_energy(mode="DMRG") # energy with DMRG
    eed = fc.gs_energy(mode="ED") # energy with exact diag.
    return edmrg,eed

nf = 10 # number of different spinless fermionic orbitals
e_dmrg,e_ed = get_energies(nf,seed=1)
print("Energy with DMRG",e_dmrg)
print("Energy with ED",e_ed)

# sweep the number of orbitals (cheap, keeps the seed fixed per point so
# each Hamiltonian realization is reproducible)
nfs = [4,6,8,10]
edmrg_sweep = []
eed_sweep = []
for nfi in nfs:
    edi,eedi = get_energies(nfi,seed=1)
    edmrg_sweep.append(edi)
    eed_sweep.append(eedi)
    print("nf =",nfi,"DMRG =",edi,"ED =",eedi)

import matplotlib.pyplot as plt
plt.plot(nfs,edmrg_sweep,marker="o",label="DMRG")
plt.plot(nfs,eed_sweep,marker="x",linestyle="--",label="ED")
plt.xlabel("Number of orbitals")
plt.ylabel("Ground-state energy")
plt.legend()
plt.show()









