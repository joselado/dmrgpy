# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')


## Bond dimension energy convergence for an S=1/2 Heisenberg chain
from dmrgpy import spinchain
import numpy as np

n= 30 # size of the chain
spins = ["S=1/2" for i in range(n)] # S=1/2 chain
sc = spinchain.Spin_Chain(spins) # create spin chain object
h = 0 # initialize Hamiltonian
for i in range(len(spins)-1):
  h = h + sc.Sx[i]*sc.Sx[i+1]
  h = h + sc.Sy[i]*sc.Sy[i+1]
  h = h + sc.Sz[i]*sc.Sz[i+1]

bds = range(3,20,2) # bond dimension
es,des = [],[] # storage of energies and fluctuations
for maxm in bds: # loop over bond dimension
  sc.set_hamiltonian(h) # create the Hamiltonian
  sc.maxm = maxm # set the bond dimension
  e = sc.gs_energy() # get the ground state energy
  wf = sc.get_gs() ; de = wf.dot(h*(h*wf)) # Energy square
  de = np.sqrt(np.abs(de-e**2)) # energy fluctuation
  es.append(e/n) # store energy
  des.append(de/n) # energy fluctuation

  print(e,de)


import matplotlib.pyplot as plt
import matplotlib
import numpy as np

matplotlib.rcParams.update({'font.size': 14})

plt.errorbar(bds,es,yerr=des,c="red",marker="o",capsize=5)
plt.ylabel("Energy")
plt.xlabel("Bond dimension")

plt.tight_layout()
plt.show()





