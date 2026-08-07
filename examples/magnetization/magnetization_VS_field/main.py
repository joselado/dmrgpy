# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import spinchain

n = 10 # number of sites in your chain
spins = ["S=1/2" for i in range(n)] # create the sites
sc = spinchain.Spin_Chain(spins) # create the chain

# now define the Hamiltonian
def geth(b):
  h = 0
  for i in range(n-1): 
      h = h + sc.Sx[i]*sc.Sx[i+1] # add exchange
      h = h + sc.Sy[i]*sc.Sy[i+1] # add exchange
      h = h + sc.Sz[i]*sc.Sz[i+1] # add exchange
  for i in range(n): h = h + b*sc.Sz[i] # add transverse field
  return h

Mz = 0
for i in range(n): Mz = Mz + sc.Sz[i]

bs = np.linspace(0.,4.0,20)
mz0s = [] # magnetization with DMRG
mz1s = [] # magnetization with ED
for b in bs:
  h = geth(b)
  sc.set_hamiltonian(h) # and initialize the Hamiltonian
  mz0 = sc.vev(Mz,mode="DMRG").real
  mz1 = sc.vev(Mz,mode="ED").real
  mz0s.append(mz0)
  mz1s.append(mz1)
  b = round(b,2)
  mz0 = round(mz0,2)
  mz1 = round(mz1,2)
  print("B=",b,"Mz DMRG = ",mz0,"Mz ED",mz1)

# plot the magnetization VS the applied field, comparing both solvers
plt.plot(bs,mz0s,marker="o",label="DMRG")
plt.scatter(bs,mz1s,label="ED",c="red")
plt.legend()
plt.xlabel("Magnetic field B")
plt.ylabel("Total magnetization $M_z$")
plt.show()











