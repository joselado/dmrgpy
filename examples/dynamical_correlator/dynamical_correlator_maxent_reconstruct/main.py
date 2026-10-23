# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Local spin spectrum S(w) of the S=1/2 Heisenberg chain, computed with
# KPM, which has a finite broadening. We reconstruct the distribution
# from it to recover sharper peaks. The spectrum is what inelastic
# neutron scattering or spin-polarized STM would measure.
import numpy as np
from dmrgpy import spinchain
from dmrgpy.reconstruct import reconstruct_distribution

n = 12 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])

# Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 20 # bond dimension
sc.kpmmaxm = 20 # bond dimension of the KPM moments
sc.nsweeps = 10 # number of DMRG sweeps

es = np.linspace(-0.5, 5.0, 600) # energies
name = (sc.Sx[0], sc.Sx[0]) # the two operators of the correlator
(x0, y0) = sc.get_dynamical_correlator(name=name, es=es, delta=0.1)
y0 = np.abs(y0) # the spectrum is positive
(x1, y1) = reconstruct_distribution(x0, y0, n=8, bnds=None)

import matplotlib.pyplot as plt
plt.plot(x0, y0, label="KPM")
plt.plot(x1, y1.real, "--", label="Reconstructed")
plt.xlabel("energy w")
plt.ylabel("S(w), first site")
plt.legend()
plt.show()
