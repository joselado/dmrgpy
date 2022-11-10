# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain
n = 6
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Chain(spins) # create the spin chain
h = 0 # generate a Heisenberg Hamiltonian
for i in range(n-1):
    h = h +sc.Sx[i]*sc.Sx[i+1]
    h = h +sc.Sy[i]*sc.Sy[i+1]
    h = h +sc.Sz[i]*sc.Sz[i+1]

h = h*h + 1j

print("Trace ED",sc.trace(h,mode="ED"))
print("Trace MPS",sc.trace(h,mode="DMRG"))

print("Inverse trace ED",sc.inverse_trace(h,mode="ED"))
print("Inverse trace MPS",sc.inverse_trace(h,mode="DMRG"))











