# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Spin spectrum S(w) of Sx_0 in the S=1/2 Heisenberg chain, computed from
# the real-time evolution of the ground state (submode="TD"). The time
# evolution uses TEBD (nearest-neighbor chains, pure-Python backend only)
# and we compare it with the KPM spectrum. Peaks are the spin
# excitations seen in inelastic neutron scattering. The peaks agree, the
# widths differ a little because the two methods broaden in different ways.
import numpy as np
from dmrgpy import spinchain

n = 5 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
sc.setup_python() # TEBD exists only in the pure-Python backend
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 30 # bond dimension

es = np.linspace(-1.0, 10.0, 1000) # energies
name = (sc.Sx[0], sc.Sx[0])
x1, y1 = sc.get_dynamical_correlator(es=es, name=name, submode="KPM", delta=0.15)
sc.tevol_method = "TEBD"
x2, y2 = sc.get_dynamical_correlator(es=es, name=name, submode="TD", delta=0.15)

import matplotlib.pyplot as plt
plt.plot(x1, y1.real, linewidth=4, label="KPM")
plt.plot(x2, np.abs(y2), "--", linewidth=2, label="TD (TEBD)")
plt.xlabel("energy w")
plt.ylabel("S(w), first site")
plt.legend()
plt.show()
