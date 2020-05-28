# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
n = 2
fc = fermionchain.Fermionic_Chain(n) # create the chain

h = fc.C[0]*fc.C[1]
h = h + h.get_dagger()
fc.set_hamiltonian(h)
print(fc.get_excited())

e0 = fc.gs_energy(mode="ED") # energy with exact diagonalization
e1 = fc.gs_energy(mode="DMRG") # energy with DMRG
print("Energy with ED",e0)
print("Energy with DMRG",e1)


