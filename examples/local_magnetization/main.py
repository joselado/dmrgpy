# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
import matplotlib.pyplot as plt
####################################
### Create the spin chain object ###
####################################
n = 10 # total number of spins
spins = [2 for i in range(n)] # list with the different spins of your system
# the spins are labeled by 2s+1, so that 2 means s=1/2, 3 meand S=1 ....
sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain object


##############################
### Create the hamiltonian ###
##############################
def fj(i,j): # function to define the exchange couplings
    if abs(i-j)==1: return 1.0  # frist neighbors
    else: return 0.0 # otherwise

bs = np.random.random((sc.ns,3))

def fb(i): # function to add local magnetic fields
    return bs[i]
    if i==0: return [0.,0.,0.4] # field only for the first atom
    else: return [0.,0.,0.] # otherwise
sc.set_exchange(fj) # add the exchange couplings
sc.set_fields(fb) # add local magnetic fields



# this parameters control the DMRG algorithm, in principle default ones are fine
#sc.maxm = 40 # maximum bond dimension
#sc.nsweeps = 12 # number of DMRG sweeps
e0 = sc.gs_energy() # compute ground state energy
print("Ground state energy is",e0)
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


