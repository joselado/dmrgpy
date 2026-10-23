# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Local spin spectrum of a 4-site Heisenberg chain from complex-time
# evolution (submode "TDZ"): the correlator is evolved in complex time
# and the real-axis result is rebuilt with a Taylor expansion in alpha0.
# We compare it with exact diagonalization (ED).
import numpy as np
from dmrgpy import spinchain

n = 4 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
sc.setup_python() # pure Python backend
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 20 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps

es = np.linspace(-0.5, 4.0, 200) # energies
delta = 0.15 # broadening of the peaks
name = (sc.Sz[0], sc.Sz[0])
(x1, y1) = sc.get_dynamical_correlator(mode="ED", submode="ED", name=name, es=es, delta=delta)
(x2, y2) = sc.get_dynamical_correlator(mode="DMRG", submode="TDZ", name=name, es=es,
        delta=delta, alpha0=0.1, n_max=4, dt=0.05)
print("Max difference TDZ-ED =", np.max(np.abs(np.array(y2) - np.array(y1))))

import matplotlib.pyplot as plt
plt.plot(x1, np.real(y1), linewidth=4, label="ED")
plt.plot(x2, np.real(y2), "--", linewidth=2, label="TDZ")
plt.xlabel("energy w")
plt.ylabel("S(w), first site")
plt.legend()
plt.show()
