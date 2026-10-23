# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Transverse field Ising chain, H = sum_i Sz_i Sz_{i+1} + Bx sum_i Sx_i.
# The fluctuation of the magnetization, <Mx^2> - <Mx>^2 with
# Mx = sum_i Sx_i / N, changes strongly at the quantum phase transition.
# It is computed with DMRG as a function of the field Bx.
import numpy as np
from dmrgpy import spinchain

N = 12 # number of sites
bxs = np.linspace(0., 1., 20) # transverse fields
fluct = []
for bx in bxs:
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(N)])
    h = 0
    for i in range(N-1): h = h + sc.Sz[i]*sc.Sz[i+1]
    for i in range(N): h = h + bx*sc.Sx[i]
    sc.set_hamiltonian(h)
    sc.maxm = 20 # bond dimension
    sc.nsweeps = 10 # number of DMRG sweeps
    wf = sc.get_gs()
    mx = sum(sc.Sx)/N # magnetization along x
    m1 = wf.dot(mx*wf).real
    m2 = wf.dot(mx*(mx*wf)).real
    fluct.append(m2 - m1**2)

import matplotlib.pyplot as plt
plt.plot(bxs, fluct, "o-")
plt.xlabel("transverse field Bx")
plt.ylabel("<Mx^2> - <Mx>^2")
plt.show()
