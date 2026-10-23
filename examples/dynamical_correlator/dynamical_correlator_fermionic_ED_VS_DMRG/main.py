# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Spectral function of a fermionic chain, A(w) = sum_n |<n|c_0|GS>|^2 delta(w - E_n + E_0)
# (c_0 acts as the second operator of the pair). It is the quantity
# measured by tunneling spectroscopy or photoemission.
# We compute it with KPM in DMRG, with KPM in ED and with the
# exact inversion (INV) in ED, and compare them.
import numpy as np
from dmrgpy import fermionchain

n = 4 # number of sites
fc = fermionchain.Fermionic_Chain(n)
np.random.seed(1)
h = 0
for i in range(n-1): # random complex hopping
    h = h + (np.random.random() + 1j*np.random.random())*fc.Cdag[i]*fc.C[i+1]
for i in range(n): # onsite energy
    h = h + np.random.random()*fc.Cdag[i]*fc.C[i]
h = h + h.get_dagger()
fc.set_hamiltonian(h)
fc.maxm = 40 # bond dimension
fc.nsweeps = 10 # number of DMRG sweeps

name = (fc.Cdag[0], fc.C[0]) # the two operators of the correlator
es = np.linspace(-3.0, 6.0, 300) # energies
delta = 0.1 # broadening of the peaks
x0, y0 = fc.get_dynamical_correlator(mode="ED", submode="INV", name=name, es=es, delta=delta)
x1, y1 = fc.get_dynamical_correlator(mode="ED", submode="KPM", name=name, es=es, delta=delta)
x2, y2 = fc.get_dynamical_correlator(mode="DMRG", submode="KPM", name=name, es=es, delta=delta)
print("Max difference DMRG-ED (KPM) =", np.max(np.abs(y2-y1)))

import matplotlib.pyplot as plt
plt.plot(x0, y0.real, linewidth=4, label="ED inversion")
plt.plot(x1, y1.real, "--", linewidth=2, label="ED KPM")
plt.plot(x2, y2.real, ":", linewidth=2, label="DMRG KPM")
plt.xlabel("energy w")
plt.ylabel("A(w)")
plt.legend()
plt.show()
