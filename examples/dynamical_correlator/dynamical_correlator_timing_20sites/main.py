# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Cost of the four DMRG methods for the dynamical correlator of the
# S=1/2 Heisenberg chain: KPM, CVM, TD and TDZ. All of them give the
# local spin spectrum S(w), here for Sz on the first site, and we measure
# the time each one takes on the same ground state.
import time
import numpy as np
from dmrgpy import spinchain

n = 12 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
h = 0
for i in range(n-1): # Heisenberg Hamiltonian
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 20 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps
sc.get_gs() # solve the ground state once

es = np.linspace(0.0, 3.0, 10) # energies
delta = 0.15 # broadening
name = (sc.Sz[0], sc.Sz[0])
submodes = ["KPM", "CVM", "TD", "TDZ"]
times = []
for submode in submodes:
    t0 = time.time()
    sc.get_dynamical_correlator(mode="DMRG", submode=submode, name=name, es=es, delta=delta)
    times.append(time.time()-t0)
    print(submode, "time (s) =", times[-1])

import matplotlib.pyplot as plt
plt.bar(submodes, times)
plt.xlabel("submode")
plt.ylabel("time (s)")
plt.show()
