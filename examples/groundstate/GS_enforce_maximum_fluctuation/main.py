# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import spinchain
from dmrgpy import fermionchain
n = 10
spins = ["S=1" for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Chain(spins) # create the spin chain
h = 0
for i in range(n-1):
    h = h +sc.Sx[i]*sc.Sx[i+1]
    h = h +sc.Sy[i]*sc.Sy[i+1]
    h = h +sc.Sz[i]*sc.Sz[i+1]

sc.maxm = 10

# Sweep the requested maximum energy fluctuation (maxde) and check that
# gs_energy() actually enforces it -- the achieved fluctuation
# (gs_energy_fluctuation()) should track the requested tolerance.
maxdes = [1e-1, 1e-3, 1e-6]
energies, fluctuations = [], []
for maxde in maxdes:
    sc.set_hamiltonian(h) # reset so each maxde is enforced from scratch
    e = sc.gs_energy(maxde=maxde) # compute the ground state energy
    de = sc.gs_energy_fluctuation()
    print("maxde",maxde,"Energy",e,"fluctuation",de)
    energies.append(e)
    fluctuations.append(de)

plt.loglog(maxdes, fluctuations, "o-", label="achieved fluctuation")
plt.loglog(maxdes, maxdes, "k--", label="requested maxde")
plt.xlabel("requested maxde")
plt.ylabel("achieved ground state energy fluctuation")
plt.title("GS_enforce_maximum_fluctuation: n=%d S=1 chain, maxm=%d" % (n, sc.maxm))
plt.legend()
plt.grid(alpha=0.3, which="both")
plt.tight_layout()
plt.show()









