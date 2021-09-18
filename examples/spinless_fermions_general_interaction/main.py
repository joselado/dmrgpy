# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')
# for example PATH = /home/jose/programs/dmrgpy/src
#PATH = PATH_TO_DMRGPY_LIBRARY

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
n = 6 # number of different spinless fermionic orbitals

# fc is an object that contains the information of the many body system
fc = fermionchain.Fermionic_Chain(n) # create the object

# create a random Hermitian hopping matrix
m = np.matrix(np.random.random((n,n)) + 1j*np.random.random((n,n)))
m = m + m.H # make it Hermitian


# add the hoppings to the fc object, the single particle part will be
# t(i,j) c^\dagger_i c_j
def t(i,j): return m[i,j] # define the function
fc.set_hoppings(t) # add the term to the Hamiltonian

# now define your many body interactions through a function V(i,j,k,l)
# the code is implemented so that the interaction added is
# V(i,j,k,l) c^\dagger_i c_j c^\dagger_k c_l
# as example we will implement here first neighbor density interaction

def vijkl(i,j,k,l):
    """Function defining the many body interaction"""
    if i==j and k==l and abs(i-k)==1: return 1.0
    else: return 0.0
fc.set_vijkl(vijkl) # add interaction term


# These parameters control the accuracy of the DMRG algorithm
fc.maxm = 30 # bond dimension, increase it for higher accuracy
fc.nsweeps = 30 # number of iterations for the MPS optimization
##################



print("GS energy with ED",fc.gs_energy(mode="ED")) # energy with exact diag
print("GS energy with DMRG",fc.gs_energy(mode="DMRG")) # energy with DMRG

# now we can compute excited state energies
#print("Excited states with ED",
#        fc.get_excited(n=4,mode="ED")) # energies with exact diag
#print("Excited states with DMRG",
#        fc.get_excited(n=4,mode="DMRG")) # energies with DMRG


####################################################

# now compute some ground state correlators
# the definition of the correlator is
# <GS|c^\dagger_i c_j |GS>
# where |GS> is the many body ground state
# give the list of pairs (ij) that you want to compute

nc = 10 # compute ten random correlators
pairs = [fc.Cdag[np.random.randint(n)]*fc.C[np.random.randint(n)] for i in range(nc)]
cs0 = np.array([fc.vev(p,mode="DMRG") for p in pairs])
cs1 = np.array([fc.vev(p,mode="ED") for p in pairs])

# now plot real and imaginary part
plt.subplot(121)
plt.scatter(range(nc),cs0.real,c="red",s=200,label="DMRG")
plt.plot(range(nc),cs1.real,c="blue",marker="o",label="ED")
plt.xlabel("correlator #")
plt.ylabel("real")
plt.legend()


plt.subplot(122)
plt.scatter(range(nc),cs0.imag,c="red",s=200,label="DMRG")
plt.plot(range(nc),cs1.imag,c="blue",marker="o",label="ED")
plt.xlabel("correlator #")
plt.ylabel("imaginary")
plt.legend()



plt.show()

# We can now compute 













