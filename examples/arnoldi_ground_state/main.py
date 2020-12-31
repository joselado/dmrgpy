# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
n = 6
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Chain(spins) # create the spin chain
h = 0
for i in range(n-1):
    h = h +np.random.random()*sc.Sx[i]*sc.Sx[i+1]
    h = h +np.random.random()*sc.Sy[i]*sc.Sy[i+1]
    h = h +np.random.random()*sc.Sz[i]*sc.Sz[i+1]

sc.set_hamiltonian(h)
#print(sc.get_excited(n=40,mode="ED")) ; exit()
from dmrgpy import mpsalgebra
e,wf = mpsalgebra.mpsarnoldi(sc,h,mode="GS")
print("Arnoldi",e.real)
print("DMRG",sc.gs_energy())
