# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Spin spectrum of a ring of three S=1/2 spins with Heisenberg exchange.
# We add the three components, S(w) = sum_a <S^a_i S^a_i>(w), and compute
# it with three different ED methods (CVM, INV and ED). They must give
# the same spectrum.
import numpy as np
from dmrgpy import spinchain
from dmrgpy.dynamicstk import spincorrelators

n = 3 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
h = 0
for i in range(n): # ring, the last site couples to the first
    j = (i+1) % n
    h = h + sc.Sx[i]*sc.Sx[j] + sc.Sy[i]*sc.Sy[j] + sc.Sz[i]*sc.Sz[j]
h = h + h.get_dagger()
sc.set_hamiltonian(h)
sc.get_gs(mode="ED") # ground state

es = np.linspace(-0.5, 6.0, 100) # energies
import matplotlib.pyplot as plt
for submode, size in [("CVM", 12), ("INV", 6), ("ED", 3)]:
    (x, y) = spincorrelators.get_full_SS_correlator(sc, mode="ED", submode=submode,
            i=0, es=es, delta=1e-2)
    plt.plot(x, y.real, "o", markersize=size, label=submode)
plt.xlabel("energy w")
plt.ylabel("S(w), first site")
plt.legend()
plt.show()
