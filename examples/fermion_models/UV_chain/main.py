# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import fermionchain

def get_energies(n):
    """DMRG and ED ground-state energies of the UV (Hubbard + nearest
    neighbor density-density) chain with n spinful sites."""
    fc = fermionchain.Spinful_Fermionic_Chain(n) # create the chain
    h = 0
    for i in range(n-1): # hopping
        h = h + fc.Cdagup[i]*fc.Cup[i+1]
        h = h + fc.Cdagdn[i]*fc.Cdn[i+1]
    for i in range(n): # Hubbard
        h = h + (fc.Nup[i]-.5)*(fc.Ndn[i]-.5)
    for i in range(n-1): # V-interaction
        h = h + (fc.Nup[i]+fc.Ndn[i])*(fc.Nup[i+1]+fc.Ndn[i+1])
    h = h + h.get_dagger()
    ##############################
    # Setup the Many Body Hamiltonian
    fc.maxm = 40
    fc.nsweeps = 40
    fc.set_hamiltonian(h) # set the hoppings
    return fc.gs_energy(mode="DMRG"),fc.gs_energy(mode="ED")

n = 4 # number of spinful fermionic sites
e_dmrg,e_ed = get_energies(n)
print(e_dmrg)
print(e_ed)

# sweep the chain length (small sizes stay cheap for this model)
ns = [2,3,4,5]
edmrg_sweep,eed_sweep = [],[]
for ni in ns:
    edi,eedi = get_energies(ni)
    edmrg_sweep.append(edi)
    eed_sweep.append(eedi)
    print("n =",ni,"DMRG =",edi,"ED =",eedi)

import matplotlib.pyplot as plt
plt.plot(ns,edmrg_sweep,marker="o",label="DMRG")
plt.plot(ns,eed_sweep,marker="x",linestyle="--",label="ED")
plt.xlabel("Number of sites")
plt.ylabel("Ground-state energy")
plt.legend()
plt.show()
