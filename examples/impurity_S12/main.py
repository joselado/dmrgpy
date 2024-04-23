# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
n = 60
spins = ["S=1/2" for i in range(n)]  
# spin 1/2 heisenberg chain
sc = spinchain.Spin_Chain(spins) # create the spin chain
h = 0
for i in range(n-1):
    h = h +sc.Sx[i]*sc.Sx[i+1]
    h = h +sc.Sy[i]*sc.Sy[i+1]
    h = h +sc.Sz[i]*sc.Sz[i+1]

B = 1.0 # local field
h = h - B*sc.Sz[0] # local perturbation

sc.set_hamiltonian(h)

mz = [sc.vev(sc.Sz[i]) for i in range(n)]

import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 22

plt.plot(range(len(mz)),mz,c="red",marker="o")
#plt.rcParams['figure.figsize'] = [20, 8]
r = np.max(np.abs(mz))*1.1
plt.ylim([-r,r])
plt.xlabel("n")
plt.ylabel("$\\langle S^z_n \\rangle$")


plt.tight_layout()
plt.show()





