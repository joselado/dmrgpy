# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Spin spectrum S(w) of the half-filled Hubbard chain at large U, compared
# with the S=1/2 Heisenberg chain. At large U the charge is frozen and the
# spins interact with the exchange J = 4t^2/U, so the two spectra coincide
# once the energy is measured in units of J. Inelastic neutron scattering
# or spin-polarized STM measure this spectrum.
import numpy as np
from dmrgpy import fermionchain, spinchain

n = 4 # number of sites
U = 10.0 # Hubbard interaction, in units of the hopping t
es = np.linspace(-1, 6, 400) # energies, in units of t

# Hubbard chain
fc = fermionchain.Spinful_Fermionic_Chain(n)
h = 0
for i in range(n-1):
    h = h + fc.Cdagup[i]*fc.Cup[i+1] + fc.Cdagdn[i]*fc.Cdn[i+1]
h = h + h.get_dagger() # hopping
for i in range(n):
    h = h + U*(fc.Nup[i]-0.5)*(fc.Ndn[i]-0.5)
fc.set_hamiltonian(h)
J = 4./U # effective exchange
(e1, s1) = fc.get_dynamical_correlator(name=(fc.Sx[1], fc.Sx[1]), es=es, mode="DMRG")

# Heisenberg chain
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
(e2, s2) = sc.get_dynamical_correlator(name=(sc.Sx[1], sc.Sx[1]), es=es/J,
                                       delta=0.1/J, mode="DMRG")

import matplotlib.pyplot as plt
plt.plot(e2, s2.real, linewidth=4, label="Heisenberg")
plt.plot(e1/J, s1.real*J, "--", linewidth=2, label="Hubbard")
plt.xlabel("energy w / J")
plt.ylabel("S(w), site 1")
plt.legend()
plt.show()
