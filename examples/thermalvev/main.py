# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
n = 2
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Chain(spins) # create the spin chain
h = 0
J,B = 0,-1

for i in range(n-1):
    h = h +J*sc.Sx[i]*sc.Sx[i+1]
    h = h +J*sc.Sy[i]*sc.Sy[i+1]
    h = h +J*sc.Sz[i]*sc.Sz[i+1]
for i in range(n): 
    h = h + B*sc.Sz[i]
sc.set_hamiltonian(h) # set the Hamiltonian

Ts = np.linspace(1e-7,4.,10) # temperatures
Szs = [sc.vev(sc.Sz[0],mode="ED",T=T) for T in Ts] # expectation values

import matplotlib.pyplot as plt

plt.plot(Ts,Szs)
plt.xlabel("Temperature")
plt.ylabel("Magnetization")
plt.show()







