# Add the root path of the dmrgpy library and uncomment the next two lines
# for example PATH = /home/jose/programs/dmrgpy/src
#PATH = PATH_TO_DMRGPY_LIBRARY
#import os ; import sys ; sys.path.append(PATH)

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
n = 12 # number of different spinless fermionic orbitals


# create a random Hermitian hopping matrix
m = np.matrix(np.random.random((n,n)) + 1j*np.random.random((n,n)))
m = m + m.H # make it Hermitian


# add the hoppings to the fc object, the single particle part will be
# t(i,j) c^\dagger_i c_j
def t(i,j): 
    if abs(i-j)==1: return 1.0
    return m[i,j] # define the function

# add now the function for hubbard
def fh(i,j):
    if abs(i-j)==1: return 1.0
    return 0.0


# create a chain and add hubbard
fc = fermionchain.Spinful_Fermionic_Hamiltonian(n//2) # create the object
fc.set_hoppings(t) # add the term to the Hamiltonian
fc.set_hubbard(fh) # add the term to the Hamiltonian
e0 = fc.gs_energy() # now get the energy



print("Energy with DMRG",fc.gs_energy(mode="DMRG"))
print("Energy with ED",fc.gs_energy(mode="ED"))
print("Magnetization with DMRG",fc.get_magnetization(mode="DMRG"))
print("Magnetization with ED",fc.get_magnetization(mode="ED"))


