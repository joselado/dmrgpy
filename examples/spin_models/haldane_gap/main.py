# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import matplotlib.pyplot as plt
from dmrgpy import spinchain
n = 18 # take n sites
spins = ["S=1" for i in range(n)] # spin 1 Heisenberg chain
sc = spinchain.Spin_Chain(spins) # create the spin chain
sc.maxm = 20
sc.nsweeps = 30
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1]
    h = h + sc.Sy[i]*sc.Sy[i+1]
    h = h + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
# compute n excited states
es = sc.get_excited(n=6) # compute excited states with DMRG
# for the Haldane chain with open boundary conditions,
# in the limit of infnite length we expect a four-fold degeneracy
# (that correspond to the dangling S=1/2 on the edges)
# so that the bulk gap is computed as the energy difference
# between 1st and the 5th state
gap = es[4] - es[0]
print("Haldane gap =",gap)
print("Energies",es-es[0])

def haldane_gap(ni,nsweeps,ne):
    """Compute the Haldane gap of an open S=1 Heisenberg chain of length ni"""
    spinsi = ["S=1" for i in range(ni)]
    sci = spinchain.Spin_Chain(spinsi)
    sci.maxm = 20
    sci.nsweeps = nsweeps
    hi = 0
    for i in range(ni-1):
        hi = hi + sci.Sx[i]*sci.Sx[i+1]
        hi = hi + sci.Sy[i]*sci.Sy[i+1]
        hi = hi + sci.Sz[i]*sci.Sz[i+1]
    sci.set_hamiltonian(hi)
    esi = sci.get_excited(n=ne)
    return esi[4] - esi[0]

# sweep the chain length: the Haldane gap should approach its known
# thermodynamic-limit value (~0.41 J) as n grows -- the canonical way to
# visualize the Haldane-gap phenomenon. Use fewer sweeps/excited states
# for the extra (smaller, cheaper-to-converge) sizes and reuse the
# already-computed, more accurate n=18 point above rather than
# recomputing it.
ns_extra = [8,14]
gaps_extra = [haldane_gap(ni,nsweeps=8,ne=5) for ni in ns_extra]
ns = ns_extra + [n]
gaps = gaps_extra + [gap]

plt.plot(ns,gaps,marker="o")
plt.xlabel("Chain length")
plt.ylabel("Haldane gap")
plt.show()
