# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Distribution of the magnetization M = (1/n) sum_i Sz_i in the ground
# state of a transverse field Ising chain, P(m) = sum_n |<n|M|GS>|^2 delta(m - m_n).
# We compute it with the kernel polynomial method (KPM): first the Chebyshev
# moments, then we rebuild the function from them with a chosen kernel.
import numpy as np
from dmrgpy import spinchain
from dmrgpy.algebra.kpm import reconstruct_chebyshev

n = 20 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])

# transverse field Ising Hamiltonian, H = -sum_i Sz_i Sz_{i+1} - 0.3 sum_i Sx_i
h = 0
for i in range(n-1):
    h = h - sc.Sz[i]*sc.Sz[i+1]
for i in range(n):
    h = h - 0.3*sc.Sx[i]
sc.set_hamiltonian(h)
sc.maxm = 20 # bond dimension of the ground state
sc.kpmmaxm = sc.maxm # bond dimension of the KPM vectors

M = (1/n)*sum(sc.Sz) # normalized magnetization
mus, shift, scale = sc.get_distribution_moments(X=M, delta=2e-2) # moments
ms = np.linspace(-0.7, 0.7, 2000) # values of the magnetization
x, y = reconstruct_chebyshev(mus, shift=shift, scale=scale, x=ms, kernel="jackson")

import matplotlib.pyplot as plt
plt.plot(x, y.real)
plt.xlabel("magnetization m")
plt.ylabel("distribution P(m)")
plt.show()
