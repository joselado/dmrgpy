# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import fermionchain
n = 6
fc = fermionchain.Fermionic_Chain(n) # create the spin chain
h = 0
for i in range(n-1):
    h = h +fc.Cdag[i]*fc.C[i+1]

h = h + h.get_dagger()

for i in range(n):
    h = h + (-1)**i*fc.Cdag[i]*fc.C[i]*0.6j # add some imaginary part


fc.set_hamiltonian(h) # set Hamiltonian



for mode in ["ED","DMRG"]:
    print("Computing using",mode,"mode")
    es = fc.get_excited(mode=mode,n=4)
    for e in es:
        print("Energies",np.round(e,2))


