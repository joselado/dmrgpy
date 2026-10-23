# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Transverse field Ising chain, H = -sum_i Sz_i Sz_{i+1} + B sum_i Sx_i.
# The ferromagnetic order, the magnetization per site <Sz>, is large
# for small B and disappears when the transverse field is increased.
# We start every DMRG run from the ferromagnetic state, so the chain
# stays in one of the two degenerate ground states.
import numpy as np
from dmrgpy import spinchain

n = 40 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
Mz = 0 # total magnetization
for i in range(n): Mz = Mz + sc.Sz[i]
sc.maxm = 20 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps
sc.set_hamiltonian(-Mz)
wf0 = sc.get_gs().copy() # ferromagnetic state, initial guess

bs = np.linspace(-1.0, 1.0, 21) # transverse fields
mzs = []
for b in bs:
    h = 0
    for i in range(n-1): h = h - sc.Sz[i]*sc.Sz[i+1]
    for i in range(n): h = h + b*sc.Sx[i]
    sc.set_hamiltonian(h)
    sc.set_initial_wf_guess(wf0)
    mzs.append(sc.vev(Mz).real/n) # magnetization per site

import matplotlib.pyplot as plt
plt.plot(bs, mzs, "o-")
plt.xlabel("transverse field B")
plt.ylabel("<Sz> per site")
plt.show()
