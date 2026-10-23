# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Full spin spectrum of a closed S=1/2 Heisenberg ring,
# S(w) = sum_n <GS|S_i.|n><n|.S_i|GS> delta(w - E_n + E_0), with S.S = Sx Sx + Sy Sy + Sz Sz.
# It adds the three spin directions, as in an unpolarized neutron experiment.
# We plot it for the first site of the ring.
from dmrgpy import spinchain

n = 5 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
h = 0
for i in range(n): # the last bond closes the ring
    j = (i+1) % n
    h = h + sc.Sx[i]*sc.Sx[j] + sc.Sy[i]*sc.Sy[j] + sc.Sz[i]*sc.Sz[j]
sc.set_hamiltonian(h)

import numpy as np
es = np.linspace(-0.2, 3, 200) # energies
(x, y) = sc.get_full_SS_correlator(es=es, delta=0.1, i=0)

import matplotlib.pyplot as plt
plt.plot(x, y.real, linewidth=3)
plt.xlabel("energy w")
plt.ylabel("S(w), site 0")
plt.show()
