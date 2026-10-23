# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Single-magnon band E(k) of an infinite XX chain in a field h above
# saturation. The ground state is the fully polarized product state, so the
# exact band is E(k) = 2h - J cos(k). We compute it with the quasiparticle
# ansatz on top of VUMPS. The band is seen in inelastic neutron scattering.
import numpy as np
from dmrgpy import infinitechain

J, h = 1.0, 3.0 # exchange and field, h > J
ic = infinitechain.Infinite_Spin_Chain(["1/2"])
ic.gs_method = "vumps" # needed for the excitations
H = -J*(ic.SxC[0]*ic.SxR[0] + ic.SyC[0]*ic.SyR[0]) - 2*h*ic.SzC[0]
ic.set_hamiltonian(H)
ic.maxm = 4 # bond dimension
ic.maxiter = 300 # maximum number of VUMPS iterations
ic.etol = 1e-12 # energy tolerance
ic.gs_energy()

ks = np.linspace(0, 2*np.pi, 13) # momenta
exact = 2*h - J*np.cos(ks)
band = np.array([ic.excitation_energies(k, n=1)[0].real for k in ks])
print("Max difference =", np.max(np.abs(band-exact)))
print("Gap =", ic.excitation_gap(), "exact =", 2*h-J)
assert np.allclose(band, exact, atol=1e-6)

import matplotlib.pyplot as plt
plt.plot(ks, exact, linewidth=4, label="exact")
plt.plot(ks, band, "o", label="VUMPS")
plt.xlabel("momentum k")
plt.ylabel("E(k)")
plt.legend()
plt.show()
