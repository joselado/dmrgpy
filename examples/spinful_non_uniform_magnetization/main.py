# Add the root path of the dmrgpy library and uncomment the next two lines
# for example PATH = /home/jose/programs/dmrgpy/src
#PATH = PATH_TO_DMRGPY_LIBRARY
#import os ; import sys ; sys.path.append(PATH)

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
from pygra import geometry
n = 40 # number of different spinless fermionic orbitals (2 times sites)
g = geometry.bichain()
g = g.supercell(n//4) # as many sites
print("Chain with ",len(g.r),"sites")
g.dimensionality = 0 # zero dimensional
h = g.get_hamiltonian(has_spin=True)
h.add_antiferromagnetism(lambda r: 1.0*(r[0]<0))
h.add_rashba(0.2)



# add the hoppings to the fc object, the single particle part will be
# t(i,j) c^\dagger_i c_j
def t(i,j): 
    return h.intra[i,j]

# create a chain and add hubbard
fc = fermionchain.Spinful_Fermionic_Hamiltonian(n//2) # create the object
fc.maxm = 20
fc.set_hoppings(t) # add the term to the Hamiltonian

if n<10: # perform the test
  print("Energy with DMRG",fc.gs_energy(mode="DMRG"))
  print("Energy with ED",fc.gs_energy(mode="ED"))
  print("Magnetization with DMRG",fc.get_magnetization(mode="DMRG"))
  print("Magnetization with ED",fc.get_magnetization(mode="ED"))


import matplotlib.pyplot as plt

mz = fc.get_magnetization()[:,2]
mz2 = h.compute_vev(name="sz").real
inds = range(len(mz))
plt.scatter(inds,mz,c="red",label="DMRG")
plt.plot(inds,mz2,c="green",label="Exact")
plt.legend()
plt.xlabel("Site")
plt.ylabel("Magnetization")
plt.show()
print(mz)



