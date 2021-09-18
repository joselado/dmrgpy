# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np # conventional numpy library
import matplotlib.pyplot as plt # library to plot the results

from dmrgpy import spinchain
ns = 6 # number of sites in the spin chain
spins = [3 for i in range(ns)] # S=1 chain
sc = spinchain.Spin_Chain(spins) # create spin chain object

# now define a custom Hamiltonian
h = 0.0 # initialize Hamiltonian
Si = [sc.Sx,sc.Sy,sc.Sz] # store the three components
for i in range(ns-1): # loop 
    for S in Si: h = h + S[i]*S[i+1]  # bilinear
    for S in Si: h = h + 1./3.*S[i]*S[i+1]*S[i]*S[i+1]  # biquadratic
sc.set_hamiltonian(h) # create the Hamiltonian
print("Energy with DMRG",sc.gs_energy(mode="DMRG"))
print("Energy with ED",sc.gs_energy(mode="ED"))










