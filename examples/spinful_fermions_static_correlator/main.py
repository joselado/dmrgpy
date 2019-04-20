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


def step(i):
    return 1.0
    return (np.tanh(g.r[i][0])+1.0)/2.

#h.add_antiferromagnetism(lambda r: step(g.get_index(r)))

def funu(i):
    """
    Function to compute the Hubbard U
    """
    return 2.0*step(i)



# add the hoppings to the fc object, the single particle part will be
# t(i,j) c^\dagger_i c_j
def t(i,j): 
    if i==j:
      return h.intra[i,j] - 2*funu(i//2) # shift in the chemical potential
    else:
      return h.intra[i,j]

# add now the function for hubbard, spinful sites
def fh(i,j):
    if i==j: return funu(i)
    return 0.0


# create a chain and add hubbard
fc = fermionchain.Spinful_Fermionic_Hamiltonian(n//2) # create the object
fc.maxm = 20
fc.kpmmaxm = fc.maxm
fc.set_hoppings(t) # add the term to the Hamiltonian
fc.set_hubbard_spinful(fh) # add the term to the Hamiltonian
fc.set_swave_pairing(lambda i: 0.01+0.4*(1-step(i))) # add the term to the Hamiltonian
e0 = fc.gs_energy() # now get the energy


if n<16: # perform the test
  print("Energy with DMRG",fc.gs_energy(mode="DMRG"))
  print("Energy with ED",fc.gs_energy(mode="ED"))
  print("Density with DMRG",fc.get_density(mode="DMRG"))
  print("Density with ED",fc.get_density(mode="ED"))



import matplotlib.pyplot as plt

pairs = [(0,i) for i in range(n//2)]

ds = fc.get_correlator(pairs=pairs,name="cdc")
ps = fc.get_correlator(pairs=pairs,name="deltadelta")
mz = fc.get_correlator(pairs=pairs,name="ZZ")
fs = fc.get_correlator(pairs=pairs,name="densitydensity") -1.0

inds = np.array(range(len(ds)))+1

ds *= inds
ps *= inds
mz *= inds
fs *= inds


def splot(n):
  plt.subplot(n)
  plt.xlabel("Site")

plt.figure()
# do all the plots
splot(221)
plt.plot(inds,ds,c="red",label="creation correlator")
plt.legend()

splot(222)
plt.plot(inds,ps,c="blue",label="pairing correlator")
plt.legend()

splot(223)
plt.plot(inds,mz,c="black",label="ZZ correlator")
plt.legend()

splot(224)
plt.plot(inds,fs,c="green",label="density correlator")
plt.legend()


plt.show()



