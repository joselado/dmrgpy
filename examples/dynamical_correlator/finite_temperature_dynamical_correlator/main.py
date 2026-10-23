# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Spin spectrum of an S=1 Heisenberg chain at zero and at finite
# temperature, S(w) = sum_nm (p_n - p_m) |<m|Sz_i|n>|^2 delta(w - E_m + E_n).
# At finite T excited states are populated, so new peaks appear at low
# energy. This is what inelastic neutron scattering sees on warming.
# The finite temperature calculation is done with exact diagonalization.
import numpy as np
from dmrgpy import spinchain

n = 6 # number of sites
sc = spinchain.Spin_Chain(["S=1" for i in range(n)])
h = 0
for i in range(n-1): # Heisenberg Hamiltonian
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)

i = 2 # site of the correlator
es = np.linspace(-1.0, 4.0, 300) # energies
(x1, y1) = sc.get_dynamical_correlator(mode="ED", i=i, j=i, name="ZZ", es=es, delta=0.05)
(x2, y2) = sc.get_dynamical_correlator(mode="ED", T=0.3, submode="ED", i=i, j=i, name="ZZ", es=es, delta=0.05)

import matplotlib.pyplot as plt
plt.plot(x1, np.abs(y1), label="T=0")
plt.plot(x2, np.abs(y2), "--", label="T=0.3")
plt.xlabel("energy w")
plt.ylabel("S(w)")
plt.legend()
plt.show()
