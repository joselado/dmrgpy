# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')
# for example PATH = /home/jose/programs/dmrgpy/src
#PATH = PATH_TO_DMRGPY_LIBRARY

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
n = 6 # number of different spinless fermionic orbitals


# create a random Hermitian hopping matrix
m = np.matrix(np.random.random((n,n)) + 1j*np.random.random((n,n)))
m = m + m.H # make it Hermitian


# add the hoppings to the fc object, the single particle part will be
# t(i,j) c^\dagger_i c_j
def t(i,j): return m[i,j] # define the function


# now define the function for generalized interaction
def vijkl(i,j,k,l):
    """Function defining the many body interaction"""
    if i==j and k==l and abs(i-k)==1: return 1.0
    else: return 0.0

# add now the function for hubbard
def fh(i,j):
    if abs(i-j)==1: return 1.0
    return 0.0


# create a chain and add hubbard
fc = fermionchain.Fermionic_Chain(n) # create the object
fc.set_hoppings(t) # add the term to the Hamiltonian
fc.set_hubbard(fh) # add the term to the Hamiltonian
e0 = fc.gs_energy() # now get the energy


# clean teh Hamiltonian and add vijkl interaction
fc = fermionchain.Fermionic_Chain(n) # create the object
fc.set_hoppings(t) # add the term to the Hamiltonian
fc.set_vijkl(vijkl) # add the vijkl interaction
e1 = fc.gs_energy() # now get the energy

print("Energy with Hubbard",e0)
print("Energy with Vijkl",e1)










