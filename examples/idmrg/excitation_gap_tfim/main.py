# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Single-magnon dispersion E(k) of the infinite transverse-field Ising chain,
# H = -sum_i X_i X_{i+1} - g sum_i Z_i, with Pauli matrices X and Z.
# For g>1 the chain is a gapped paramagnet, and the exact result of the
# free fermions is E(k) = 2 sqrt(1 + g^2 - 2 g cos(k)), with gap 2(g-1).
# We compute the ground state with VUMPS and the excitations on top of it.
import numpy as np
from dmrgpy import infinitechain

g = 2.5 # transverse field
ic = infinitechain.Infinite_Spin_Chain(["1/2"])
ic.gs_method = "vumps"
# Sx = X/2 and Sz = Z/2, so that H = -4 Sx Sx - 2 g Sz
ic.set_hamiltonian(-4.0*ic.SxC[0]*ic.SxR[0] - 2.0*g*ic.SzC[0])
ic.maxm = 2 # bond dimension
ic.etol = 1e-12 # energy tolerance
ic.vumps_nrestarts = 6
print("Energy density", ic.gs_energy())

ks = np.linspace(0, np.pi, 13) # momenta
exact = 2*np.sqrt(1 + g**2 - 2*g*np.cos(ks))
got = np.array([ic.excitation_energies(k, n=1)[0].real for k in ks])
print("Max difference VUMPS-exact =", np.max(np.abs(got-exact)))
assert np.allclose(got, exact, atol=2e-3)

import matplotlib.pyplot as plt
plt.plot(ks, exact, linewidth=4, label="exact")
plt.plot(ks, got, "o--", linewidth=2, label="VUMPS")
plt.xlabel("momentum k")
plt.ylabel("E(k)")
plt.legend()
plt.show()
