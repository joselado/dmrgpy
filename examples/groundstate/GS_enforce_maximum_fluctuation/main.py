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
# gs_energy() actually enforces it -- the achieved fluctuation per site
# (gs_energy_fluctuation()/n, since maxde is a fluctuation per site) should
# sit at or below the requested tolerance. The last point sits just above
# it: each retry doubles maxm and runs two sweeps, the loop stops after five
# retries, and by maxm=160 the energy is exact to every printed digit
# (-12.894560, ED -12.894560) with 1.8e-6 per site left.
maxdes = [1e-1, 1e-3, 1e-6]
energies, fluctuations = [], []
for maxde in maxdes:
    sc.set_hamiltonian(h) # reset so each maxde is enforced from scratch
    e = sc.gs_energy(maxde=maxde) # compute the ground state energy
    de = sc.gs_energy_fluctuation()/n # per site, the unit maxde is in
    print("maxde",maxde,"Energy",e,"fluctuation per site",de)
    energies.append(e)
    fluctuations.append(de)

plt.loglog(maxdes, fluctuations, "o-", label="achieved fluctuation per site")
plt.loglog(maxdes, maxdes, "k--", label="requested maxde")
plt.xlabel("requested maxde")
plt.ylabel("achieved energy fluctuation per site")
plt.title("GS_enforce_maximum_fluctuation: n=%d S=1 chain, maxm=%d" % (n, sc.maxm))
plt.legend()
plt.grid(alpha=0.3, which="both")
plt.tight_layout()
plt.show()









