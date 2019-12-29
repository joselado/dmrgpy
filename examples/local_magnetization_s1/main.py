# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
import matplotlib.pyplot as plt
####################################
### Create the spin chain object ###
####################################
n = 11 # total number of spins
spins = [3 for i in range(n)] # list with the different spins of your system
# the spins are labeled by 2s+1, so that 2 means s=1/2, 3 means S=1 ....
sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain object
##############################
### Create the hamiltonian ###
##############################
def fj(i,j): # function to define the exchange couplings
    if abs(i-j)==1: return 1.0  # first neighbors
    if i==j: # anisotropic term, given as a matrix
        # the anisotropy can be introduced by providing a matrix m_nm
        # for the exchange, that translates into couplings of the
        # form m_nm S^n_i S^m_j, with n,m the components and
        # i j the indexes of the sites (for onsite anisotropy i=j)
        d = np.zeros((3,3)) # zero matrix
        d[0,0] = -0.3 # add easy axis anisotrpy in the x direction
        return d
    else: return 0.0 # otherwise
def fb(i): # function to add local magnetic fields
    if i==n//2: return [0.,0.,0.4] # field only for the central atom
    else: return [0.,0.,0.] # otherwise
sc.set_exchange(fj) # add the exchange couplings
sc.set_fields(fb) # add local magnetic fields
# this parameters control the DMRG algorithm, in principle default ones are fine
sc.maxm = 40 # maximum bond dimension
sc.nsweeps = 24 # number of DMRG sweeps
e0 = sc.gs_energy() # compute ground state energy
print("Ground state energy is",e0)
mx,my,mz = sc.get_magnetization() # get expectation value of the magnetization
plt.plot(range(len(mz)),mz,marker="o")
plt.xlabel("Site #")
plt.ylabel("Expectation value of Sz")
plt.show()
# uncomment this to check the code (only for small systems!!)
#e1 = sc.gs_energy(mode="ED") # ground state energy using exact diagonalization
#print("Energy using exact diagonalization",e1)


