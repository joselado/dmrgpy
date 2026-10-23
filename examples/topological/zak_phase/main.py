# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Zak phase of a dimerized Heisenberg chain, J_i = 1 + delta (-1)^i. A twist
# exp(2 pi i k) in the last bond gives a ground state |k>. The Zak phase is
# the phase of the product <k_1|k_2><k_2|k_3>...<k_N|k_1>. It is 0 or pi
# and it changes when the dimerization changes sign (topological transition).
import numpy as np
from dmrgpy import spinchain

n = 4 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])

def zak(delta):
    wfs = []
    for k in np.linspace(0., 1., 20, endpoint=False): # twist angle
        h = 0
        for i in range(n-1):
            J = 1. + delta*(-1)**i # exchange
            h = h + J*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
        z = np.exp(2j*np.pi*k) # twist
        Spn = sc.Sx[n-1] + 1j*sc.Sy[n-1] # S+ in the last site
        Sp0 = sc.Sx[0] + 1j*sc.Sy[0] # S+ in the first site
        J = 1. + delta*(-1)**(n-1)
        h = h + J*(z*Spn*Sp0.get_dagger()/2 + np.conj(z)*Sp0*Spn.get_dagger()/2
                   + sc.Sz[0]*sc.Sz[n-1]) # last bond with the twist
        sc.set_hamiltonian(h)
        wfs.append(sc.get_gs(mode="DMRG"))
    phi = 1.
    for i in range(len(wfs)): # product of overlaps, closing the loop
        phi = phi*wfs[i].dot(wfs[(i+1) % len(wfs)])
    return (np.angle(phi)/np.pi + 0.5) % 2 - 0.5 # phase / pi, 0 or 1

deltas = np.linspace(-0.6, 0.6, 13) # dimerization
phis = [zak(d) for d in deltas]

import matplotlib.pyplot as plt
plt.plot(deltas, phis, "o-")
plt.xlabel("Dimerization delta")
plt.ylabel("Zak phase / pi")
plt.show()
