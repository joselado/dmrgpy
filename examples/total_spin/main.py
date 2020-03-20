# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np # conventional numpy library
import matplotlib.pyplot as plt # library to plot the results

from dmrgpy import spinchain
ns = 6 # number of sites in the spin chain
spins = [3 for i in range(ns-1)]+[2] # S=1 chain
sc = spinchain.Spin_Chain(spins) # create spin chain object

# now define a custom Hamiltonian
h = 0.0 # initialize Hamiltonian
Si = [sc.Sx,sc.Sy,sc.Sz] # store the three components
for i in range(ns-1): # loop 
    for S in Si: h = h + S[i]*S[i+1]  # bilinear
sc.set_hamiltonian(h) # create the Hamiltonian

St = 0
sxt = 0
syt = 0
szt = 0
for s in sc.Sx: sxt = sxt + s
for s in sc.Sy: syt = syt + s
for s in sc.Sz: szt = szt + s

print("Total spin DMRG",sc.vev(sxt*sxt+syt*syt+szt*szt,mode="DMRG"))
print("Total spin ED",sc.vev(sxt*sxt+syt*syt+szt*szt,mode="ED"))

