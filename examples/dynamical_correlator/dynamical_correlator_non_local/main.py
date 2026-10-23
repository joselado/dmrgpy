# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Non-local spin correlator of the S=1/2 Heisenberg chain,
# S(r,w) = sum_n <GS|Sz_r|n><n|Sz_0|GS> delta(w - E_n + E_0).
# It tells how a spin excitation created at site 0 shows up at site r.
# Its Fourier transform in r is the S(q,w) seen in neutron scattering.
import numpy as np
from dmrgpy import spinchain

n = 6 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 20 # bond dimension
sc.kpmmaxm = 20 # bond dimension of the KPM
sc.nsweeps = 10 # number of DMRG sweeps

es = np.linspace(-0.5, 3.0, 200) # energies
zs = []
for r in range(n): # correlator between site 0 and site r
    name = (sc.Sz[r], sc.Sz[0])
    (e, s) = sc.get_dynamical_correlator(mode="DMRG", name=name, es=es, delta=0.05)
    zs.append(s.real)
zs = np.array(zs).T # energy along the first axis

import matplotlib.pyplot as plt
plt.contourf(range(n), es, zs, 100, cmap="rainbow")
plt.xlabel("distance r")
plt.ylabel("energy w")
plt.colorbar(label="S(r,w)")
plt.show()
