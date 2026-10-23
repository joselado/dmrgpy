# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Local spin spectrum S(w) of the S=1/2 Heisenberg chain, with the
# correction vector computed away from the real axis (CVMimag). The
# spectrum is then continued to the real axis with a Pade approximant.
# This smooths the peaks, so we compare with exact diagonalization (ED).
# S(w) is what inelastic neutron scattering would measure.
import numpy as np
from dmrgpy import spinchain

n = 6 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])

# Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 40 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps

es = np.linspace(-0.5, 3.0, 200) # energies
name = (sc.Sx[0], sc.Sx[0]) # the two operators of the correlator
(x1, y1) = sc.get_dynamical_correlator(mode="ED", submode="ED", name=name, es=es, delta=0.05)
(x2, y2) = sc.get_dynamical_correlator(submode="CVMimag", name=name, es=es, delta=0.05)
w1, w2 = np.trapezoid(y1.real, x1), np.trapezoid(y2.real, x2)
print("Spectral weight ED =", w1, " CVMimag =", w2)

import matplotlib.pyplot as plt
plt.plot(x1, y1.real, linewidth=4, label="ED")
plt.plot(x2, y2.real, "--", linewidth=2, label="CVMimag")
plt.xlabel("energy w")
plt.ylabel("S(w), first site")
plt.legend()
plt.show()
