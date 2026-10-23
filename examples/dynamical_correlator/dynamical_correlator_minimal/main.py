# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Minimal dynamical correlator, S(w) = sum_n |<n|Sz_0|GS>|^2 delta(w-E_n+E_0),
# of the S=1/2 Heisenberg chain. We compute it with DMRG (KPM) and with
# exact diagonalization (ED), and compare both.
import numpy as np
from dmrgpy import spinchain

n = 6 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 20 # bond dimension
sc.kpmmaxm = 20 # bond dimension in the KPM

name = (sc.Sz[0], sc.Sz[0]) # the two operators of the correlator
es = np.linspace(-0.5, 4.5, 500) # energies
delta = 0.05 # broadening
(x1, y1) = sc.get_dynamical_correlator(mode="DMRG", name=name, es=es, delta=delta)
(x2, y2) = sc.get_dynamical_correlator(mode="ED", submode="ED", name=name, es=es, delta=delta)

import matplotlib.pyplot as plt
plt.plot(x2, y2.real, linewidth=4, label="ED")
plt.plot(x1, y1.real, "--", linewidth=2, label="DMRG")
plt.xlabel("energy w")
plt.ylabel("S(w), first site")
plt.legend()
plt.show()
