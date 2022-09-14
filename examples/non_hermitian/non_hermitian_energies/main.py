# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain
n = 10
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Chain(spins) # create the spin chain
h = 0
for i in range(n-1):
    h = h +sc.Sx[i]*sc.Sx[i+1]
    h = h +sc.Sy[i]*sc.Sy[i+1]
    h = h +sc.Sz[i]*sc.Sz[i+1]

for i in range(n):
    h = h + (-1)**i*sc.Sz[i]*0.3j # add some imaginary part

#h = h +0.2j

sc.set_hamiltonian(h) # set Hamiltonian

for mode in ["DMRG"]:
    print("Computing using",mode,"mode")
    es = sc.get_excited(mode=mode,n=4,verbose=1)
    for e in es:
        print("Energies",np.round(e,2))


