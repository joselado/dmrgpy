# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Zak phase of a dimerized spinful fermion chain (SSH chain), with hoppings
# 1 +- delta. A flux k in the boundary hopping of spin up gives a ground state |k>,
# and the Zak phase is the phase of the product of overlaps <k_i|k_{i+1}>
# around the loop. It is 0 or pi, and it jumps at delta=0.
import numpy as np
from dmrgpy import fermionchain

n = 4 # number of sites
fc = fermionchain.Spinful_Fermionic_Chain(n)

def zak(delta):
    ks = np.linspace(0., 1., 20, endpoint=False) # fluxes
    wfs = []
    for k in ks:
        h = 0
        for i in range(n):
            t = 1. + delta*(-1)**i # dimerized hopping
            j = (i+1) % n
            z = np.exp(2j*np.pi*k) if j == 0 else 1. # flux in the last bond
            h = h + t*z*fc.Cdagup[i]*fc.Cup[j] + t*fc.Cdagdn[i]*fc.Cdn[j]
        fc.set_hamiltonian(h + h.get_dagger())
        wfs.append(fc.get_gs(mode="ED"))
    phi = wfs[-1].dot(wfs[0]) # close the loop
    for i in range(len(wfs)-1):
        phi = phi*wfs[i].dot(wfs[i+1])
    return np.angle(phi)/np.pi

deltas = np.linspace(-0.6, 0.6, 13)
phis = [zak(d) for d in deltas]

import matplotlib.pyplot as plt
plt.plot(deltas, phis, "o-")
plt.xlabel("dimerization delta")
plt.ylabel("Zak phase / pi")
plt.show()
