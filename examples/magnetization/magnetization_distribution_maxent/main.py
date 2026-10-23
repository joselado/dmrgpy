# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Distribution of the staggered magnetization M = sum_i (-1)^i Sz_i / (2n)
# in the ground state of an S=1 chain with single-ion anisotropy.
# The distribution P(m) gives the probability of each value of M in a
# projective measurement. We compute it with the kernel polynomial method
# (KPM) and with the maximum entropy method, and compare them.
import numpy as np
from dmrgpy import spinchain
from dmrgpy import distribution

n = 12 # number of sites
sc = spinchain.Spin_Chain(["S=1" for i in range(n)])
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    h = h - 1.1*sc.Sz[i]*sc.Sz[i] # single-ion anisotropy
sc.set_hamiltonian(h)
sc.maxm = 20 # bond dimension
sc.kpmmaxm = sc.maxm # bond dimension of the KPM vectors
sc.nsweeps = 10 # number of DMRG sweeps

M = sum([sc.Sz[i]*(-1)**i for i in range(n)])/(2*n)
x, y = sc.get_distribution(X=M, delta=5e-2) # KPM
x1, y1 = distribution.get_distribution_maxent(sc, X=M, n=6) # max entropy
print("Integral KPM =", np.trapezoid(y.real, dx=x[1]-x[0]))
print("Integral max ent =", np.trapezoid(y1.real, dx=x1[1]-x1[0]))

import matplotlib.pyplot as plt
plt.plot(x, y.real, label="KPM")
plt.plot(x1, y1.real, "--", label="Max ent")
plt.xlabel("magnetization")
plt.ylabel("distribution")
plt.legend()
plt.show()
