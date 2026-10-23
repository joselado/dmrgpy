# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Phase transition of the S=1 Heisenberg ring with single-ion anisotropy,
# H = sum_i S_i . S_{i+1} + D sum_i Sz_i^2. For small D the ground state is
# a Haldane state, and for large D it is a trivial state with Sz=0 in all
# sites. We follow the transition with the fluctuation of the
# magnetization M = sum_i Sx_i / n in the ground state, <M^2> - <M>^2.
import numpy as np
from dmrgpy import spinchain

n = 6 # number of sites
sc = spinchain.Spin_Chain(["S=1" for i in range(n)])
M = sum(sc.Sx)/n # magnetization

def fluctuation(D):
    h = 0
    for i in range(n): # periodic ring
        j = (i+1) % n
        h = h + sc.Sx[i]*sc.Sx[j] + sc.Sy[i]*sc.Sy[j] + sc.Sz[i]*sc.Sz[j]
        h = h + D*sc.Sz[i]*sc.Sz[i]
    sc.set_hamiltonian(h)
    sc.maxm = 40 # bond dimension
    wf = sc.get_gs() # ground state
    return (wf.dot(M*(M*wf)) - wf.dot(M*wf)**2).real

Ds = np.linspace(0., 2., 16) # anisotropies
fs = [fluctuation(D) for D in Ds]

import matplotlib.pyplot as plt
plt.plot(Ds, fs, "o-")
plt.xlabel("single-ion anisotropy D")
plt.ylabel("<M^2> - <M>^2")
plt.show()
