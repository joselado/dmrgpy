# Add the root path of the dmrgpy library and uncomment the next two lines
# for example PATH = /home/jose/programs/dmrgpy/src
#PATH = PATH_TO_DMRGPY_LIBRARY
#import os ; import sys ; sys.path.append(PATH)

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
from pygra import geometry
n = 16 # number of different spinless fermionic orbitals (2 times sites)
g = geometry.bichain()
g = g.supercell(n//4) # as many sites
print("Chain with ",len(g.r),"sites")
g.dimensionality = 0 # zero dimensional
h = g.get_hamiltonian(has_spin=True)


def step(i):
    return (np.tanh(g.r[i][0])+1.0)/2.

h.add_antiferromagnetism(lambda r: step(g.get_index(r)))

def funu(i):
    """
    Function to compute the Hubbard U
    """
    return 0.0*step(i)



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
fc.set_hubbard(fh) # add the term to the Hamiltonian
fc.set_swave_pairing(lambda i: 0.4*(1-step(i))) # add the term to the Hamiltonian
e0 = fc.gs_energy() # now get the energy


if n<16: # perform the test
  print("Energy with DMRG",fc.gs_energy(mode="DMRG"))
  print("Energy with ED",fc.gs_energy(mode="ED"))
  print("Density with DMRG",fc.get_density(mode="DMRG"))
  print("Density with ED",fc.get_density(mode="ED"))



import matplotlib.pyplot as plt

ds = fc.get_density()
ps = fc.get_pairing()
mz = fc.get_magnetization()[:,2]
fs = fc.get_density_fluctuation()
ps = ps/(ps[0]/np.abs(ps[0])) # factor out the phase
ps = ps.real
inds = range(len(ds))


# compute the dynamical correlator
if n<20:
  delta = 0.1
  print("Computing dynamical correlator")
  (es,dcs) = fc.get_dynamical_correlator(name="cdc",i=0,j=0,delta=delta)
  (es,dcb) = fc.get_dynamical_correlator(name="cdc",i=n//4,j=n//4,delta=delta)
  (es,dca) = fc.get_dynamical_correlator(name="cdc",i=n//2-1,j=n//2-1,delta=delta)
  plt.figure()
  plt.plot(es,dcs,label="SC",c="red")
  plt.plot(es,dcb,label="interface",c="blue")
  plt.plot(es,dca,label="AF",c="green")
  plt.legend()
  print("Finished dynamical correlator")


def splot(n):
  plt.subplot(n)
  plt.xlabel("Site")

plt.figure()
# do all the plots
splot(221)
plt.scatter(inds,ds,c="red",label="density DMRG")
plt.legend()

splot(222)
plt.scatter(inds,ps,c="blue",label="pairing DMRG")
plt.legend()

splot(223)
plt.scatter(inds,mz,c="black",label="magnetization DMRG")
plt.legend()

splot(224)
plt.plot(inds,fs,c="green",label="density fluctuation DMRG")
plt.legend()


plt.show()



