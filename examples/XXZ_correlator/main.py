# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
n = 20
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Chain(spins) # create the spin chain
h = 0 # 
D = 0.5
for i in range(n-1):
    h = h +sc.Sx[i]*sc.Sx[i+1]
    h = h +sc.Sy[i]*sc.Sy[i+1]
    h = h +D*sc.Sz[i]*sc.Sz[i+1]

sc.set_hamiltonian(h)
sxi = [sc.vev(sc.Sx[0]*sc.Sx[i]) for i in range(n)]
syi = [sc.vev(sc.Sy[0]*sc.Sy[i]) for i in range(n)]
szi = [sc.vev(sc.Sz[0]*sc.Sz[i]) for i in range(n)]

import matplotlib.pyplot as plt

plt.plot(range(len(sxi)),sxi,label="XX",marker="o",markersize=10)
plt.plot(range(len(syi)),syi,label="YY",marker="o",markersize=6)
plt.plot(range(len(szi)),szi,label="ZZ",marker="o",markersize=6)
plt.xlabel("Distance") ; plt.ylabel("Correlator")
plt.legend()
plt.show()







