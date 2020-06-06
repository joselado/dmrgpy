# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
import matplotlib.pyplot as plt
####################################
### Create the spin chain object ###
####################################
n = 6 # total number of spins
spins = ["S=1/2" for i in range(n)] # list with the different spins 
# the spins are labeled by 2s+1, so that 2 means s=1/2, 3 meand S=1 ....
sc = spinchain.Spin_Chain(spins) # create the spin chain object
sc.itensor_version = "julia"

##############################
### Create the hamiltonian ###
##############################
h = 0 
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1]
    h = h + sc.Sy[i]*sc.Sy[i+1]
    h = h + sc.Sz[i]*sc.Sz[i+1]

h = h + 0.1*sc.Sz[0]
sc.set_hamiltonian(h)
#sc.itensor_version = "julia"
# this parameters control the DMRG algorithm, in principle default ones are fine
#sc.maxm = 40 # maximum bond dimension
#sc.nsweeps = 12 # number of DMRG sweeps
e0 = sc.gs_energy() # compute ground state energy
e1 = sc.gs_energy(mode="ED") # compute ground state energy
print("Energy with DMRG is",e0)
print("Energy with ED is",e1)
mx,my,mz = sc.get_magnetization() # use DMRG
mx2,my2,mz2 = sc.get_magnetization(mode="ED") # use ED

# plot the resutls
plt.plot(range(len(mz)),mz,label="DMRG",c="red")
plt.scatter(range(len(mz2)),mz2,label="ED",c="blue")
plt.legend()
plt.xlabel("Site #")
plt.ylabel("Expectation value of Sz")
plt.show()
# uncomment this to check the code (only for small systems!!)
#e1 = sc.gs_energy(mode="ED") # ground state energy using exact diagonalization
#print("Energy using exact diagonalization",e1)


