# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Distribution of the magnetization M = sum_i Sz_i / n in the ground state
# of a Heisenberg chain with a transverse field. The ground state is not
# an eigenstate of M, so P(m) = sum_k |<k|GS>|^2 delta(m - m_k) has weight
# at several values. We compute it with the kernel polynomial method (KPM)
# and compare two kernels, which differ in how they broaden the peaks.
import numpy as np
from dmrgpy import spinchain

n = 6 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
for i in range(n): h = h - 0.6*sc.Sx[i] # transverse field
sc.set_hamiltonian(h)
sc.maxm = 20 # bond dimension
sc.kpmmaxm = sc.maxm
sc.nsweeps = 10 # number of DMRG sweeps

M = sum(sc.Sz)/n # operator whose distribution we want
(x1, y1) = sc.get_distribution(X=M, delta=5e-2, kernel="lorentz")
(x2, y2) = sc.get_distribution(X=M, delta=5e-2, kernel="jackson")
print("Weight, Lorentz =", np.trapezoid(y1.real, dx=x1[1]-x1[0]))

import matplotlib.pyplot as plt
plt.plot(x1, y1.real, label="Lorentz")
plt.plot(x2, y2.real, label="Jackson")
plt.xlabel("magnetization")
plt.ylabel("distribution")
plt.legend()
plt.show()
