# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Local spin spectrum S(w) of the S=1/2 Heisenberg chain, computed from the
# time evolution (submode="TD") and with the kernel polynomial method (KPM).
# The time signal is damped, which gives the broadening delta, and with
# predict=True it is extended in time by linear prediction.
# We compare both choices of predict with KPM. In this example they overlap,
# and the time evolution has longer tails than KPM.
import numpy as np
from dmrgpy import spinchain

n = 6 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 30 # bond dimension

es = np.linspace(-0.5, 3.0, 400) # energies
name = (sc.Sz[0], sc.Sz[0]) # the two operators of the correlator
delta = 0.1 # broadening of the peaks
(x0, y0) = sc.get_dynamical_correlator(es=es, name=name, submode="KPM", delta=delta)
(x1, y1) = sc.get_dynamical_correlator(es=es, name=name, submode="TD", delta=delta, predict=False)
(x2, y2) = sc.get_dynamical_correlator(es=es, name=name, submode="TD", delta=delta, predict=True)

import matplotlib.pyplot as plt
plt.plot(x0, np.abs(y0), linewidth=4, label="KPM")
plt.plot(x1, np.abs(y1), "--", linewidth=2, label="TD, predict=False")
plt.plot(x2, np.abs(y2), "--", linewidth=2, label="TD, predict=True")
plt.xlabel("energy w")
plt.ylabel("S(w), first site")
plt.legend()
plt.show()
