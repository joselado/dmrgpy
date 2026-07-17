# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Dynamical correlator computed via real-time evolution of an MPS
# (submode="TD", dynamics.py -> timedependent.dynamical_correlator() ->
# evolution_dmrg_DC()), which now defaults to the TDVP engine
# (mpscpp3/TDVP/) instead of the legacy MPO-Taylor method. Compared here
# against the (unrelated, KPM-based) submode="KPM" method as a sanity
# check that the two independent methods agree on the main spectral
# features.
import numpy as np
from dmrgpy import spinchain

n = 5
spins = ["S=1/2" for i in range(n)]
sc = spinchain.Spin_Chain(spins) # create the chain
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1]
    h = h + sc.Sy[i]*sc.Sy[i+1]
    h = h + sc.Sz[i]*sc.Sz[i+1]

sc.set_hamiltonian(h)

sc.maxm = 10
es = np.linspace(-1.0,10.0,2000)
name = (sc.Sx[0],sc.Sx[0]) # correlator to compute

print("Time-evolution method for the TD submode:",sc.tevol_method) # TDVP

import time
print("Starting")
t0 = time.time()
(x1,y1) = sc.get_dynamical_correlator(es=es,name=name,submode="KPM")
t1 = time.time()
(x2,y2) = sc.get_dynamical_correlator(es=es,name=name,submode="TD")
t2 = time.time()
print("Time in TD (TDVP-backed)",t2-t1)
print("Time in KPM",t1-t0)

import matplotlib.pyplot as plt
plt.plot(x1,y1.real,label="KPM")
plt.plot(x2,np.abs(y2),label="TD (TDVP)")
plt.legend()
plt.xlabel("frequency")
plt.ylabel("Dynamical correlator")
plt.show()
