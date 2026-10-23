# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Local spin spectrum S(w) = sum_n |<n|Sz_0|GS>|^2 delta(w - E_n + E_0)
# of an S=1/2 chain, computed with DMRG and the kernel polynomial method.
# The kernel decides how each delta function is broadened. The Jackson
# kernel is the default, and the Lorentz kernel gives a different line shape
# for the same number of moments.
import numpy as np
from dmrgpy import spinchain

n = 8 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])

# anisotropic Heisenberg Hamiltonian, the xx coupling is 0.7
h = 0
for i in range(n-1):
    h = h + 0.7*sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 20 # bond dimension
sc.kpmmaxm = 20 # bond dimension used in KPM
sc.nsweeps = 10 # number of DMRG sweeps

name = (sc.Sz[3], sc.Sz[3]) # the two operators of the correlator
es = np.linspace(-0.5, 4.5, 1000) # energies
delta = 0.1 # broadening of the peaks
(x1, y1) = sc.get_dynamical_correlator(name=name, es=es, delta=delta)
(x2, y2) = sc.get_dynamical_correlator(name=name, es=es, delta=delta,
        kernel="lorentz")

import matplotlib.pyplot as plt
plt.plot(x1, y1.real, linewidth=3, label="Jackson kernel")
plt.plot(x2, y2.real, "--", linewidth=2, label="Lorentz kernel")
plt.xlabel("energy w")
plt.ylabel("S(w), site 3")
plt.legend()
plt.show()
