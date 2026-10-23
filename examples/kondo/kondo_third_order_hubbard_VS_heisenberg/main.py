# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Third-order Kondo spectrum dI/dV (Ternes, arXiv:1505.04430) on the first
# site of a half-filled Hubbard chain, against the Heisenberg chain it
# maps onto at large U, with J = 4t^2/U. The tip couples to the local
# spin, and STM measures the spectrum. The Hubbard chain adds charge
# fluctuations, so the curves approach the Heisenberg one as U grows.
# We use exact diagonalization, and energies are in units of J.
import numpy as np
from dmrgpy import spinchain, fermionchain

n = 6 # number of sites
kw = dict(site=0, Jrho_s=-0.05, T=0., kB=1., order=3)
x = np.linspace(-3., 3., 240) # eV/J

def spectrum(chain, J):
    _, d = chain.get_kondo_spectrum(x*J, omega0=10*J, Gamma0=1e-2*J, **kw)
    return d/(2*np.pi)

sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
h = 0
for i in range(n-1): # Heisenberg chain with J=1
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
heis = spectrum(sc, 1.)

import matplotlib.pyplot as plt
plt.plot(x, heis, "k", linewidth=4, label="Heisenberg")
for U in [8., 32.]: # Hubbard chain, t=1
    fc = fermionchain.Spinful_Fermionic_Chain_Native(n, itensor_version="python")
    h = 0
    for i in range(n-1): # hopping
        h = h - fc.Cdagup[i]*fc.Cup[i+1] - fc.Cdagdn[i]*fc.Cdn[i+1]
    h = h + h.get_dagger()
    for i in range(n): # Hubbard interaction
        h = h + U*(fc.Nup[i]-0.5)*(fc.Ndn[i]-0.5)
    fc.set_hamiltonian(h)
    fc.set_conserved_sector(Nf=n) # half filling
    plt.plot(x, spectrum(fc, 4./U), "--", label="Hubbard U/t=%g" % U)
plt.xlabel("eV / J")
plt.ylabel("dI/dV")
plt.legend()
plt.show()
