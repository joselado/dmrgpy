# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Potential term of the third-order STM tunneling spectrum of a spin
# chain in a magnetic field (Ternes, New J. Phys. 17 063016), where
# the tip couples to the spin with exchange J and also with a potential U.
# It is the part of dI/dV proportional to U, and it is computed
# with DMRG (dynamical correlator) and with exact diagonalization (ED).
import numpy as np
from dmrgpy import spinchain
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
from dmrgpy.kondospectrumtk.conductance import third_order_potential_dIdV
from dmrgpy.kondospectrumtk.potentialdc import third_order_potential_dIdV_dc

muB = 5.7883818066e-5 # Bohr magneton in eV/T
sc = spinchain.Spin_Chain(["1/2", "1/2", "1/2"], itensor_version=3)
sc.maxm = 30 # bond dimension
sc.nsweeps = 15 # number of DMRG sweeps
h = 2.0*muB*10.0*sc.Sz[0] # Zeeman field of 10 T
for i in range(2):
    h = h + 0.01*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
sc.set_hamiltonian(h)
sc.get_gs()
ks = KondoSpectrum(sc, site=0, T=0.0)

eVs = np.linspace(-2e-3, 2e-3, 41) # bias voltages
es = np.linspace(-20e-3, 20e-3, 5328) # energies of the spin spectrum
Jrho, U = 0.05, 0.2 # exchange and potential couplings
d1 = third_order_potential_dIdV(ks, eVs, Jrho, U, T0=1.0)
d2 = third_order_potential_dIdV_dc(sc, 0, eVs, Jrho, U, T0=1.0,
                                   mode="DMRG", submode="KPM", delta=2e-5, es=es)
print("Max difference =", np.max(np.abs(d2-d1)))

import matplotlib.pyplot as plt
plt.plot(eVs, d1, linewidth=4, label="ED")
plt.plot(eVs, d2, "--", linewidth=2, label="DMRG")
plt.xlabel("eV")
plt.ylabel("potential term of dI/dV")
plt.legend()
plt.show()
