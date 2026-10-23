# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Third-order Kondo tunneling spectrum dI/dV at T=0 of a magnetic impurity
# (site 0) weakly coupled to two spins, as measured by STM. We compute it
# with exact diagonalization (ED) and with DMRG plus time evolution, which
# must agree. Energies are in eV, the impurity is split by a magnetic field.
import numpy as np
from dmrgpy import spinchain

sc = spinchain.Spin_Chain(["S=1/2" for i in range(3)], itensor_version=3)
sc.maxm = 30 # bond dimension
sc.nsweeps = 15 # number of DMRG sweeps
h = 2.0*5.7883818066e-5*10.0*sc.Sz[0] # Zeeman term of the impurity
for i in range(2): # weak Heisenberg coupling
    h = h + 0.01*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
sc.set_hamiltonian(h)

eVs = np.linspace(-2e-3, 2e-3, 41) # bias voltages
kw = dict(site=0, Jrho_s=0.05, U=0.0, T=0.0, order=3, omega0=2e-3, Gamma0=2e-3)
(_, ded) = sc.get_kondo_spectrum(eVs, mode="ED", **kw)
(_, ddmrg) = sc.get_kondo_spectrum(eVs, mode="DMRG", submode="KPM", delta=2e-5,
        es=np.linspace(-3e-3, 3e-3, 800), dt2=150.0, n_t2_half=20,
        dtau=150.0, n_tau_half=20, **kw)
print("Max difference =", np.max(np.abs(ddmrg-ded)))
assert np.max(np.abs(ddmrg-ded)) < 0.02*np.max(np.abs(ded))

import matplotlib.pyplot as plt
plt.plot(eVs, ded, linewidth=4, label="ED")
plt.plot(eVs, ddmrg, "--", linewidth=2, label="DMRG")
plt.xlabel("bias eV")
plt.ylabel("dI/dV")
plt.legend()
plt.show()
