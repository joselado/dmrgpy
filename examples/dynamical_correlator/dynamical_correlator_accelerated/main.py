# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Local spin spectrum S(w) = sum_n |<n|Sx_0|GS>|^2 delta(w - E_n + E_0)
# of the S=1/2 Heisenberg chain, computed with DMRG and the kernel
# polynomial method (KPM). The option kpm_accelerate=True makes KPM
# faster, and the spectrum should be the same with and without it.
import numpy as np
from dmrgpy import spinchain

n = 6 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])

# Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
h = 0
for i in range(n-1):
    h = h + sc.SS(i, i+1)
sc.set_hamiltonian(h)
sc.maxm = 20 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps

name = (sc.Sx[0], sc.Sx[0]) # the two operators of the correlator
sc.kpm_accelerate = True # fast KPM
(x0, y0) = sc.get_dynamical_correlator(mode="DMRG", name=name)
sc.kpm_accelerate = False # slow KPM
(x1, y1) = sc.get_dynamical_correlator(mode="DMRG", name=name)
print("Max difference =", np.max(np.abs(y0 - y1)))

import matplotlib.pyplot as plt
plt.plot(x1, y1.real, linewidth=4, label="Without acceleration")
plt.plot(x0, y0.real, "--", linewidth=2, label="With acceleration")
plt.xlabel("energy w")
plt.ylabel("S(w), first site")
plt.xlim([-0.5, 4.5])
plt.legend()
plt.show()
