# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import fermionchain
n = 6
fc = fermionchain.Fermionic_Chain(n) # create the spin chain
h = 0
for i in range(n-1):
    h = h +fc.Cdag[i]*fc.C[i+1]

h = h + h.get_dagger()

for i in range(n):
    h = h + (-1)**i*fc.Cdag[i]*fc.C[i]*0.6j # add some imaginary part


fc.set_hamiltonian(h) # set Hamiltonian



results = {}
for mode in ["ED","DMRG"]:
    print("Computing using",mode,"mode")
    es = fc.get_excited(mode=mode,n=4)
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


