# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
n = 16
fc = fermionchain.Fermionic_Chain(n) # create the chain
h = 0
for i in range(n-1):
    h = h + fc.C[i]*fc.Cdag[i+1]
h = h + h.get_dagger()
fc.set_hamiltonian(h)
e0 = fc.gs_energy(mode="ED") # energy with exact diagonalization
e1 = fc.gs_energy(mode="DMRG") # energy with DMRG
print("Energy with ED",e0)
print("Energy with DMRG",e1)


