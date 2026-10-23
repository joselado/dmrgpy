# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Local spin spectrum of the S=1/2 Heisenberg chain,
# S(w) = sum_n |<n|Sx_0|GS>|^2 delta(w - E_n + E_0).
# We obtain it with two methods, the kernel polynomial method (KPM)
# and the maximum entropy method (maxent), and compare them.
# The peaks are the spin excitations seen in inelastic neutron scattering.
import numpy as np
from dmrgpy import spinchain

n = 8 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 20 # bond dimension
sc.kpmmaxm = 20 # bond dimension of the KPM
sc.nsweeps = 10 # number of DMRG sweeps

es = np.linspace(-0.5, 4.5, 300) # energies
name = (sc.Sx[0], sc.Sx[0]) # the two operators of the correlator
(x1, y1) = sc.get_dynamical_correlator(submode="KPM", name=name, es=es, delta=0.1)
(x2, y2) = sc.get_dynamical_correlator(submode="maxent", name=name, es=es, delta=0.05)

import matplotlib.pyplot as plt
plt.plot(x1, np.abs(y1), linewidth=3, label="KPM")
plt.plot(x2, np.abs(y2), "--", linewidth=2, label="Max entropy")
plt.xlabel("energy w")
plt.ylabel("S(w), first site")
plt.legend()
plt.show()
