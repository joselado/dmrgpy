# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain

n = 40 # number of sites in your chain
spins = ["S=1/2" for i in range(n)] # create the sites
sc = spinchain.Spin_Chain(spins) # create the chain

# now define the Hamiltonian
def geth(b):
  h = 0
  for i in range(n-1): 
      h = h - sc.Sz[i]*sc.Sz[i+1] # add exchange
  for i in range(n): h = h + b*sc.Sx[i] # add transverse field
  return h

# define the total spin in the z-direction
Mz = 0
for i in range(n): Mz = Mz + sc.Sz[i]

# get a fully ferromagnetic wavefunction (to be used as initial guess later on)
sc.set_hamiltonian(-Mz) ; wffe = sc.get_gs().copy() # ferromagnetic wavefunction

bs = np.linspace(-1.0,1.0,40) # list of magnetic fields
mzs = []
# loop over magnetic fields
for b in bs:
    h = geth(b) # get the Hamiltonian
    sc.set_hamiltonian(h) # and initialize the Hamiltonian
    sc.set_initial_wf_guess(wffe) # setup this wavefunction as the initial guess
    mz = sc.vev(Mz).real/n # compute the magnetization per site
    print("B=",b,"Mz = ",mz)
    mzs.append(mz) # store in an array

np.savetxt("MZ_VS_B.OUT",np.array([bs,mzs]).T) # save in a file



import matplotlib.pyplot as plt


plt.plot(bs,mzs,marker="o",c="red")
plt.xlabel("B")
plt.ylabel("Sz")
plt.show()

