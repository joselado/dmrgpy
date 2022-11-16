# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

from dmrgpy import spinchain
import numpy as np
n = 10 ; spins = ["S=1/2" for i in range(n)] # spins in each site
sc = spinchain.Spin_Chain(spins) # create spin chain object
hxy,hz = 0,0 # initialize Hamiltonian
for i in range(len(spins)-1):
  hxy = hxy + sc.Sx[i]*sc.Sx[i+1]
  hxy = hxy + sc.Sy[i]*sc.Sy[i+1]
  hz = hz + sc.Sz[i]*sc.Sz[i+1]
jzs = np.linspace(0.,3.,20) # Jz exchange couplings
es = [] # storage for energies
for jz in jzs: # loop over Jz
  sc.set_hamiltonian(hxy + jz*hz) # create the Hamiltonian
  es.append(sc.gs_energy()) # compute and store energy


import matplotlib.pyplot as plt
import matplotlib
import numpy as np

matplotlib.rcParams.update({'font.size': 14})

plt.plot(jzs,es,marker="o",c="blue")
plt.ylabel("Energy")
plt.xlabel("$J_z$")
#plt.xlim([-0.2,4])

plt.tight_layout()
plt.show()





