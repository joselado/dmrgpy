# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Sum rule of the KPM dynamical correlator. For A = Sz on a spin-1/2
# site, A*A = 1/4, so the integral of S(w) is exactly 1/4, for any chain.
# We check how well the KPM spectrum fulfills it as the bond dimension
# of the moments (kpmmaxm) grows.
import numpy as np
from dmrgpy import spinchain

n = 8 # number of sites
es = np.linspace(-9, 9, 1500) # wide window, to hold all the weight
kpmmaxms = [5, 10, 20, 40] # bond dimensions of the KPM moments
errors = []
for kpmmaxm in kpmmaxms:
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = 30 # bond dimension
    sc.nsweeps = 10 # number of DMRG sweeps
    sc.kpmmaxm = kpmmaxm
    name = (sc.Sz[0], sc.Sz[0]) # the two operators of the correlator
    (x, y) = sc.get_dynamical_correlator(mode="DMRG", name=name, es=es, delta=0.1)
    weight = np.trapezoid(np.real(y), x) # integrated weight
    errors.append(abs(weight - 0.25))
    print("kpmmaxm =", kpmmaxm, " weight =", weight)

import matplotlib.pyplot as plt
plt.plot(kpmmaxms, errors, "o-")
plt.yscale("log")
plt.xlabel("kpmmaxm")
plt.ylabel("|integral of S(w) - 1/4|")
plt.show()
