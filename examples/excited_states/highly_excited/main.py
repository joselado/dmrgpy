# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain
n = 4
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Chain(spins) # create the spin chain
h = 0
for i in range(n-1):
    h = h +np.random.random()*sc.Sx[i]*sc.Sx[i+1]
    h = h +np.random.random()*sc.Sy[i]*sc.Sy[i+1]
    h = h +np.random.random()*sc.Sz[i]*sc.Sz[i+1]

sc.set_hamiltonian(h)
#print(sc.get_excited(n=40,mode="ED")) ; exit()
from dmrgpy import mpsalgebra
e,wf = mpsalgebra.mpsarnoldi(sc,h,e=.4,mode="ShiftInv")
print("Energy",e)

# sweep the shift-invert target energy to see which excited state is
# recovered for different targets (cheap, same small n=4 chain)
targets = [0.1,0.2,0.3,0.4,0.5,0.6]
found = []
for et in targets:
    ei,wfi = mpsalgebra.mpsarnoldi(sc,h,e=et,mode="ShiftInv")
    ei = np.array(ei).flatten()[0]
    found.append(ei)
    print("Target",et,"-> recovered energy",ei)

import matplotlib.pyplot as plt
plt.plot(targets,found,marker="o",label="Recovered energy")
plt.plot(targets,targets,linestyle="--",color="gray",label="Target energy")
plt.xlabel("Target energy (shift-invert)")
plt.ylabel("Recovered eigenenergy")
plt.legend()
plt.show()









