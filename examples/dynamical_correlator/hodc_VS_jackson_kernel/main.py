# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Local spin spectrum S(w) of the Heisenberg chain from the kernel
# polynomial method (KPM). The same Chebyshev moments can be summed with
# the usual Jackson kernel or with the high-order delta-Chebyshev (HODC)
# kernel, which is less smooth but sharper. We compare both with ED.
import numpy as np
from dmrgpy import spinchain

n = 8 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])

# Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 40 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps

es = np.linspace(-0.5, 3.0, 300) # energies
delta = 0.05 # broadening of the peaks
name = (sc.Sz[0], sc.Sz[0]) # the two operators of the correlator
(e0, s0) = sc.get_dynamical_correlator(name=name, es=es, delta=delta, mode="ED")
(e1, s1) = sc.get_dynamical_correlator(name=name, es=es, delta=delta, mode="DMRG", kernel="jackson")
(e2, s2) = sc.get_dynamical_correlator(name=name, es=es, delta=delta, mode="DMRG", kernel="hodc")

import matplotlib.pyplot as plt
plt.plot(e0, s0.real, "k", linewidth=4, label="ED")
plt.plot(e1, s1.real, "--", linewidth=2, label="KPM, Jackson")
plt.plot(e2, s2.real, "--", linewidth=2, label="KPM, HODC")
plt.xlabel("energy w")
plt.ylabel("S(w), first site")
plt.legend()
plt.show()
