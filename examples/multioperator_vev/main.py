# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np # conventional numpy library
from dmrgpy import spinchain # library dealing with DMRG for spin chains
import matplotlib.pyplot as plt # library to plot the results

n = 10 # total number of spins
spins = [2 for i in range(n)] # list with the different spins of your system
# the spins are labeled by 2s+1, so that 2 means s=1/2, 3 meand S=1 ....
sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain object
h = 0 
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1]
    h = h + sc.Sy[i]*sc.Sz[i+1]
    h = h + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h) # add the exchange couplings

op = sc.Sx[0]
for i in range(1,n): op = op*sc.Sx[i]

print("VEV with DMRG",sc.vev(op,mode="DMRG").real)
print("VEV with ED",sc.vev(op,mode="ED").real)
