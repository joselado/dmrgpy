# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Local spin spectrum of the S=1/2 Heisenberg chain with DMRG,
# S(w) = sum_n |<n|Sz_0|GS>|^2 delta(w - E_n + E_0).
# It is computed with the kernel polynomial method (KPM), and the
# delta functions are broadened by delta. Each peak is a spin excitation,
# which inelastic neutron scattering or spin-polarized STM can measure.
import numpy as np
from dmrgpy import spinchain

n = 6 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])

# Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 20 # bond dimension
sc.kpmmaxm = 20 # bond dimension used in KPM
sc.nsweeps = 10 # number of DMRG sweeps

name = (sc.Sz[0], sc.Sz[0]) # the two operators of the correlator
es = np.linspace(-0.5, 4, 200) # energies
delta = 5e-2 # broadening of the peaks
(x, y) = sc.get_dynamical_correlator(name=name, es=es, delta=delta)

import matplotlib.pyplot as plt
plt.plot(x, y.real, label="DMRG")
plt.xlabel("energy w")
plt.ylabel("S(w), first site")
plt.legend()
plt.show()
