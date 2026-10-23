# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# STM dI/dV of a single S=1/2 spin in a magnetic field B, from third-order
# perturbation theory (Ternes, New J. Phys. 17 063016, Fig. 3b).
# At B=0 the exchange with the electrons gives a Kondo peak at zero bias.
# A field splits it, and steps appear at eV = +-g muB B (spin-flip energy).
# The peak height at B=0 is 1.13 in units of e^2 T0^2/h in the paper.
import numpy as np
from dmrgpy import spinchain

g, muB = 2.0, 5.7883818066e-5 # g factor and Bohr magneton in eV/T
eVs = np.linspace(-4e-3, 4e-3, 161) # bias voltages in eV

import matplotlib.pyplot as plt
for B in [0.0, 4.0, 10.0]: # magnetic fields in Tesla
    sc = spinchain.Spin_Chain(["1/2"])
    sc.set_hamiltonian(g*muB*B*sc.Sz[0]) # Zeeman term
    (eV, d) = sc.get_kondo_spectrum(eVs, site=0, Jrho_s=-0.05, T=1.0,
            order=3, omega0=20e-3) # exchange Jrho_s, temperature T in K
    plt.plot(eV*1e3, d/(2*np.pi), label="B=%g T" % B)
    if B == 0.0: print("Peak at B=0 =", np.max(d)/(2*np.pi), "(paper 1.13)")
plt.xlabel("bias eV (meV)")
plt.ylabel("dI/dV (e^2 T0^2/h)")
plt.legend()
plt.show()
