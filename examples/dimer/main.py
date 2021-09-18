# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
n = 2
spins = ["S=1" for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Chain(spins) # create the spin chain
h = 0
h = h +sc.Sx[0]*sc.Sx[1]
h = h +sc.Sy[0]*sc.Sy[1]
h = h +sc.Sz[0]*sc.Sz[1]
h = -h
c0 = -0.1
S111 = sc.Sx[0] + sc.Sy[0] + sc.Sz[0]
L111 = sc.Sx[1] + sc.Sy[1] + sc.Sz[1]
h = h + c0*S111*S111



sc.set_hamiltonian(h)
print(sc.get_excited(n=9,mode="ED"))









