# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

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
mz_map = [] # magnetization per site, for every field (for the plot below)
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
    mz_map.append(mzs) # store this field's row for the heatmap
fo.close()

import matplotlib.pyplot as plt

mz_map = np.array(mz_map) # shape (len(bs), n)
plt.imshow(mz_map,aspect="auto",origin="lower",
        extent=[0,n-1,bs[0],bs[-1]],cmap="bwr")
plt.colorbar(label="Sz")
plt.xlabel("Site")
plt.ylabel("B")
plt.show()









