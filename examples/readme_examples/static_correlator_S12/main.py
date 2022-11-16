# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

from dmrgpy import spinchain
n = 30
spins = ["S=1/2" for i in range(n)] # S=1 in each site
sc = spinchain.Spin_Chain(spins) # create spin chain object
h = 0 # initialize Hamiltonian
for i in range(len(spins)-1):
  h = h + sc.Sx[i]*sc.Sx[i+1]
  h = h + sc.Sy[i]*sc.Sy[i+1]
  h = h + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h) # create the Hamiltonian
cs = [sc.vev(sc.Sz[0]*sc.Sz[i]).real for i in range(n)]

import matplotlib.pyplot as plt
import matplotlib
import numpy as np

matplotlib.rcParams.update({'font.size': 14})

plt.plot(np.array(range(len(cs))),cs,marker="o",c="red")
plt.ylabel("$\langle S^z_0 S^z_{n} \\rangle$")
plt.xlabel("$n$")

plt.tight_layout()
plt.show()





