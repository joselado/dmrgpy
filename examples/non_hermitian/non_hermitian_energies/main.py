# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain
n = 6
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Chain(spins) # create the spin chain
h = 0
for i in range(n-1):
    h = h +sc.Sx[i]*sc.Sx[i+1]
    h = h +sc.Sy[i]*sc.Sy[i+1]
    h = h +sc.Sz[i]*sc.Sz[i+1]

for i in range(n):
    h = h + (-1)**i*sc.Sz[i]*0.6j # add some imaginary part


sc.set_hamiltonian(h) # set Hamiltonian



results = {}
for mode in ["ED","DMRG"]:
    print("Computing using",mode,"mode")
    es = sc.get_excited(mode=mode,n=4)
    for e in es:
        print("Energies",np.round(e,2))
    results[mode] = np.array(es)

# plot the excited-state energies of each mode in the complex plane
import matplotlib.pyplot as plt
colors = {"ED":"red","DMRG":"blue"}
markers = {"ED":"x","DMRG":"o"}
for mode,es in results.items():
    plt.scatter(es.real,es.imag,c=colors[mode],marker=markers[mode],
                s=80,label=mode)
plt.xlabel("Re(E)")
plt.ylabel("Im(E)")
plt.legend()
plt.show()


