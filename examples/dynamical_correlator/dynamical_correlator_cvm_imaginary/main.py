# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Local spin spectrum S(w) of the S=2 Heisenberg chain.
# We compute it with two DMRG methods, the kernel polynomial method (KPM)
# and the correction vector method at imaginary frequency (CVMimag).
# Both must give the same spectrum, as inelastic neutron scattering would measure.
import numpy as np
from dmrgpy import spinchain

n = 8 # number of sites
sc = spinchain.Spin_Chain([5 for i in range(n)]) # 5 states per site, S=2

# Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 10 # bond dimension
sc.kpmmaxm = 10 # KPM bond dimension
sc.nsweeps = 10 # number of DMRG sweeps

es = np.linspace(-0.5, 5., 300) # energies
name = (sc.Sx[0], sc.Sx[0]) # the two operators of the correlator
(x1, y1) = sc.get_dynamical_correlator(submode="KPM", name=name, es=es, delta=5e-2)
(x2, y2) = sc.get_dynamical_correlator(submode="CVMimag", name=name, es=es, delta=5e-2)

import matplotlib.pyplot as plt
plt.plot(x1, np.abs(y1), linewidth=4, label="KPM")
plt.plot(x2, np.abs(y2), "--", linewidth=2, label="CVMimag")
plt.xlabel("energy w")
plt.ylabel("S(w), first site")
plt.legend()
plt.show()
