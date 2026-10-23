# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Single-magnon dispersion E(k) of the infinite transverse-field Ising chain,
# H = -sum sigma^x_i sigma^x_{i+1} - g sum sigma^z_i, in the gapped phase g>1.
# The exact result is E(k) = 2 sqrt(1 + g^2 - 2 g cos(k)), the energy that
# inelastic neutron scattering would measure. We compute it with VUMPS and
# the excitation ansatz, with the python and the v3 backends.
import numpy as np
from dmrgpy import infinitechain

g = 2.5 # transverse field
ks = np.linspace(0, np.pi, 13) # momenta
exact = 2*np.sqrt(1 + g**2 - 2*g*np.cos(ks))

def dispersion(version):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=version)
    ic.gs_method = "vumps"
    ic.set_hamiltonian(-4.0*ic.SxC[0]*ic.SxR[0] - 2.0*g*ic.SzC[0])
    ic.maxm = 3 # bond dimension
    ic.etol = 1e-10
    ic.vumps_nrestarts = 6
    ic.gs_energy()
    return np.array([ic.excitation_energies(k, n=1)[0].real for k in ks])

e_py = dispersion("python")
e_v3 = dispersion(3)
print("Max difference python-v3 =", np.max(np.abs(e_py-e_v3)))
assert np.max(np.abs(e_py-e_v3)) < 1e-5

import matplotlib.pyplot as plt
plt.plot(ks, exact, "k-", linewidth=4, label="exact")
plt.plot(ks, e_py, "o", markersize=10, label="python")
plt.plot(ks, e_v3, "x", markersize=10, label="v3")
plt.xlabel("momentum k")
plt.ylabel("E(k)")
plt.legend()
plt.show()
