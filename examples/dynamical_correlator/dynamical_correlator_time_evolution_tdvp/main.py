# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Local spin spectrum S(w) = sum_n |<n|Sx_0|GS>|^2 delta(w - E_n + E_0) of
# the S=1/2 Heisenberg chain, computed in two independent ways. With
# submode="TD" we evolve the state in real time with TDVP and Fourier
# transform. With submode="KPM" we expand the spectrum in Chebyshev
# polynomials. Both must give the same peaks.
import numpy as np
from dmrgpy import spinchain

n = 5 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
# Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 10 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps

es = np.linspace(-1.0, 10.0, 2000) # energies
name = (sc.Sx[0], sc.Sx[0]) # the two operators of the correlator
(x1, y1) = sc.get_dynamical_correlator(es=es, name=name, submode="KPM")
(x2, y2) = sc.get_dynamical_correlator(es=es, name=name, submode="TD")

import matplotlib.pyplot as plt
plt.plot(x1, y1.real, linewidth=4, label="KPM")
plt.plot(x2, np.abs(y2), "--", linewidth=2, label="TD (TDVP)")
plt.xlabel("energy w")
plt.ylabel("S(w), first site")
plt.legend()
plt.show()
