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

bs = np.linspace(0.0,1.0,40) # list of magnetic fields
mzs = []
# loop over magnetic fields
fo = open("M_VS_SITE_VS_B.OUT","w")
for b in bs:
    h = geth(b) # get the Hamiltonian
    sc.set_hamiltonian(h) # and initialize the Hamiltonian
    sc.set_initial_wf_guess(wffe) # setup this wavefunction as the initial guess
    mzs = [sc.vev(op).real for op in sc.Sz] # compute the magnetization per site
    for i in range(n):
        fo.write(str(i)+ " ")
        fo.write(str(b)+ " ")
        fo.write(str(mzs[i])+ "\n")
    print(mzs)
fo.close()









