# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np # conventional numpy library
from dmrgpy import spinchain # library dealing with DMRG for spin chains
from dmrgpy import thermal # library dealing with DMRG for spin chains
import matplotlib.pyplot as plt # library to plot the results
####################################
### Create the spin chain object ###
####################################
n = 4 # total number of spins
spins = ["S=1/2" for i in range(n)] # spins of your system
#sc = spinchain.Spin_Chain(spins) # create the spin chain object
sc = thermal.Thermal_Spin_Chain(spins) # create the spin chain object
##############################
### Create the hamiltonian ###
##############################
h = 0 
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1]
    h = h + sc.Sy[i]*sc.Sy[i+1]
    h = h + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h) # set the Hamiltonian
sc.T = 0.01
sc.mode = "ED"

wf = sc.get_gs()
print(wf.dot(h*wf).real)









