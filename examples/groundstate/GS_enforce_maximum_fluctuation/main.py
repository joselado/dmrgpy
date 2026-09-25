# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# gs_energy(maxde=...) doubles maxm until the energy fluctuation PER SITE,
# ||(H-<H>)|psi>||/n, drops below maxde. gs_energy_fluctuation() returns the
# total, so it is divided by n here before it is compared with maxde; the
# plot is the per-site fluctuation reached against the per-site request.
import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import spinchain
n = 10
spins = ["S=1" for i in range(n)] # spin 1 Heisenberg chain
sc = spinchain.Spin_Chain(spins) # create the spin chain
h = 0
for i in range(n-1):
    h = h +sc.Sx[i]*sc.Sx[i+1]
    h = h +sc.Sy[i]*sc.Sy[i+1]
    h = h +sc.Sz[i]*sc.Sz[i+1]

# the sweep schedule, pinned: maxm is the starting bond dimension the maxde
# loop doubles from
sc.maxm, sc.nsweeps = 10, 15
sc.noise, sc.cutoff = 1e-7, 1e-12

# The achieved fluctuation per site should sit at or below the requested
# tolerance. The last point sits just above it: each retry doubles maxm and
# runs two sweeps, the loop stops after five retries, and by maxm=160 the
# energy is exact to every printed digit (-12.894560, ED -12.894560) with
# 1.8e-6 per site left.
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
plt.loglog(maxdes, maxdes, "k--", label="requested maxde (per site)")
plt.xlabel("requested maxde")
plt.ylabel(r"energy fluctuation per site $\|(H-\langle H\rangle)|\psi\rangle\|/n$")
plt.title("maxde per site: n=%d S=1 chain, starting maxm=%d" % (n, sc.maxm))
plt.legend()
plt.grid(alpha=0.3, which="both")
plt.tight_layout()
plt.savefig("GS_enforce_maximum_fluctuation.png",dpi=150)
plt.show()
