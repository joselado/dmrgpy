# Add the root path of the dmrgpy library and uncomment the next two lines
# for example PATH = /home/jose/programs/dmrgpy/src
#PATH = PATH_TO_DMRGPY_LIBRARY
#import os ; import sys ; sys.path.append(PATH)

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
n = 8 # number of different spinless fermionic orbitals


# create a random Hermitian hopping matrix
m = np.matrix(np.random.random((n,n)) + 1j*np.random.random((n,n)))
m = m + m.H # make it Hermitian


# add the hoppings to the fc object, the single particle part will be
# t(i,j) c^\dagger_i c_j
def t(i,j): 
    if abs(i//2-j//2)==1 and i%2==j%2: return 1.0
    else: return 0.0
    return m[i,j] # define the function

# add now the function for hubbard, spinful sites
def fh(i,j):
    if i==j: return 1.0
    return 0.0


# create a chain and add hubbard
fc = fermionchain.Spinful_Fermionic_Chain(n//2) # create the object
fc.set_hoppings(t) # add the term to the Hamiltonian
fc.set_hubbard(fh) # add the term to the Hamiltonian
fc.set_swave_pairing(lambda i: 0.4) # add the term to the Hamiltonian
#e0 = fc.gs_energy() # now get the energy



print("Energy with DMRG",fc.gs_energy(mode="DMRG"))
print("Energy with ED",fc.gs_energy(mode="ED"))
mdm = fc.get_density(mode="DMRG")
med = fc.get_density(mode="ED")
print("Density with DMRG",mdm)
print("Density with ED",med)

print("Difference in density",np.max(np.abs((mdm-med))))


