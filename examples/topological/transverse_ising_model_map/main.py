# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Magnetization Sz_i of the transverse field Ising chain, H = -sum Sz_i
# Sz_i+1 + B sum Sx_i, for several fields B. At B=0 the chain is a
# ferromagnet, and the transverse field reduces the order, mostly at the
# edges. The ferromagnetic state is the initial guess of each DMRG run.
import numpy as np
from dmrgpy import spinchain

n = 40 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
sc.maxm = 20 # bond dimension
sc.nsweeps = 6 # number of DMRG sweeps
sc.set_hamiltonian(-sum(sc.Sz))
wf0 = sc.get_gs().copy() # ferromagnetic state

bs = np.linspace(0.0, 1.0, 20) # transverse fields
mz = []
for b in bs:
    h = 0
    for i in range(n-1): h = h - sc.Sz[i]*sc.Sz[i+1]
    for i in range(n): h = h + b*sc.Sx[i]
    sc.set_hamiltonian(h)
    sc.set_initial_wf_guess(wf0)
    mz.append([sc.vev(op).real for op in sc.Sz])

import matplotlib.pyplot as plt
plt.imshow(np.array(mz), aspect="auto", origin="lower",
        extent=[0, n-1, bs[0], bs[-1]], cmap="bwr")
plt.colorbar(label="<Sz_i>")
plt.xlabel("site i")
plt.ylabel("transverse field B")
plt.show()
