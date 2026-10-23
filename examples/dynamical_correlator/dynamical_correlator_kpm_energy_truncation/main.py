# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Local spin spectrum of a 4-site Heisenberg chain with the kernel
# polynomial method (KPM). A narrow KPM energy window gives more
# resolution, but the Chebyshev recursion can blow up. Energy truncation
# (kpm_energy_truncate) removes the high-energy part and keeps it stable.
# The lowest peak must sit at the gap of the chain, we compare with ED.
import numpy as np
from dmrgpy import spinchain

n = 4 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
sc.setup_python() # energy truncation is only in this backend
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 20 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps

es = np.linspace(0.3, 1.0, 71) # energies
name = (sc.Sz[0], sc.Sz[0])
(x1, y1) = sc.get_dynamical_correlator(mode="ED", submode="INV", name=name, es=es, delta=0.05)
sc.kpm_scale = 0.65 # narrow KPM window
sc.kpm_energy_truncate = True # remove the high-energy components
sc.kpm_truncate_nsweeps = 10
(x2, y2) = sc.get_dynamical_correlator(mode="DMRG", submode="KPM", name=name, es=es, delta=0.05)
print("Peak ED =", x1[np.argmax(y1.real)], " Peak KPM =", x2[np.argmax(y2.real)])

import matplotlib.pyplot as plt
plt.plot(x1, y1.real, linewidth=4, label="ED")
plt.plot(x2, y2.real, "--", linewidth=2, label="KPM, energy truncation")
plt.xlabel("energy w")
plt.ylabel("S(w), first site")
plt.legend()
plt.show()
