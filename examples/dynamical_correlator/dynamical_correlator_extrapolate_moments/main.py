# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Local spin spectrum S(w) of an anisotropic Heisenberg chain with the
# kernel polynomial method (KPM). A narrow peak needs many moments. With
# extrapolation we get a narrow broadening from the moments of a broad one.
import numpy as np
from dmrgpy import spinchain

n = 10 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
h = 0
for i in range(n-1):
    h = h + 0.5*sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 20 # bond dimension
sc.kpmmaxm = 20 # bond dimension in the KPM

name = (sc.Sz[1], sc.Sz[1]) # the two operators of the correlator
es = np.linspace(-0.5, 6, 1000) # energies
delta = 0.2 # broadening
factor = 2. # the extrapolation doubles the number of moments

sc.kpm_extrapolate = False
(x1, y1) = sc.get_dynamical_correlator(name=name, es=es, delta=delta)
(x2, y2) = sc.get_dynamical_correlator(name=name, es=es, delta=2*delta)
sc.kpm_extrapolate = True
sc.kpm_extrapolate_factor = factor
sc.kpm_extrapolate_mode = "power"
(x3, y3) = sc.get_dynamical_correlator(name=name, es=es, delta=delta)

import matplotlib.pyplot as plt
plt.plot(x1, y1.real, linewidth=4, label="all moments")
plt.plot(x2, y2.real, ":", linewidth=2, label="half the moments")
plt.plot(x3, y3.real, "--", linewidth=2, label="extrapolated")
plt.xlabel("energy w")
plt.ylabel("S(w), site 1")
plt.legend()
plt.show()
