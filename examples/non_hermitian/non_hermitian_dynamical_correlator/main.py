# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import fermionchain

L = 3
n = 4*L
fc = fermionchain.Fermionic_Chain(n) # create the fermion chain
h = 0
for i in range(n-1):
    h = h + fc.Cdag[i]*fc.C[i+1]
h = h + h.get_dagger()
eta = 0.5
eta0 = 0.

for i in range(L):
    h = h + 1j*(eta0 + eta)*fc.Cdag[4*i]*fc.C[4*i]
    h = h + 1j*(eta0 - eta)*fc.Cdag[4*i+1]*fc.C[4*i+1]
    h = h + 1j*(eta0 - eta)*fc.Cdag[4*i+2]*fc.C[4*i+2]
    h = h + 1j*(eta0 + eta)*fc.Cdag[4*i+3]*fc.C[4*i+3]


fc.set_hamiltonian(h)

fc.get_gs(mode="ED")
fc.get_gs(mode="ED")

out = []
energies = np.linspace(-0.1,3.0,200)
for i in range(n):
    name = fc.C[i],fc.Cdag[i]
    (es,ds) = fc.get_dynamical_correlator(name=name,mode="ED",
            es=energies,delta=1e-1)
    out.append([es,ds])



import matplotlib.pyplot as plt
plt.subplot(1,2,1)
ie = 0
for o in out:
    (es,ds) = o[0],o[1]
    plt.plot(es,ds,label=str(ie))
    ie += 1

plt.legend()
cmap = [list(reversed(ds)) for (es,ds) in out]

plt.subplot(1,2,2)
plt.imshow(np.array(cmap).T,aspect="auto") ; plt.axis("off")


plt.show()
