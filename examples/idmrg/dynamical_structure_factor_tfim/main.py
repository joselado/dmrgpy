# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Dynamical structure factor S(k,w) of the infinite transverse-field Ising
# chain H = -sx sx - g sz, in the gapped phase g>1. Sx creates one magnon, so
# S(k,w) is a delta peak on the band E(k) = 2 sqrt(1+g^2-2g cos k), with a
# weight that we get from the tangent-space excitations (VUMPS). This is
# what inelastic neutron scattering measures.
import numpy as np
from dmrgpy import infinitechain

g = 2.5 # transverse field
ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
# Sx and Sz have eigenvalues 1/2, so these are the Pauli matrices
ic.set_hamiltonian(-4.0*ic.SxC[0]*ic.SxR[0] - 2.0*g*ic.SzC[0])
ic.gs_method = "vumps" # needed for the excitations
ic.maxm = 4 # bond dimension
ic.maxiter = 800 # maximum number of VUMPS iterations
ic.vumps_nrestarts = 6
ic.gs_energy()

ks = np.linspace(-np.pi, np.pi, 61) # momenta
es = np.linspace(1.5, 8.5, 400) # energies
ks, es, S = ic.dynamical_structure_factor("Sx", ks=ks, energies=es,
        delta=0.08, n=1) # n=1, only the one-magnon branch

band = [ic.spectral_weights("Sx", k, n=1)[0][0] for k in ks]
exact = 2.0*np.sqrt(1.0 + g**2 - 2.0*g*np.cos(ks))
print("Max difference =", np.max(np.abs(np.array(band)-exact)))
assert np.max(np.abs(np.array(band)-exact)) < 1e-4

import matplotlib.pyplot as plt
plt.pcolormesh(ks, es, S.T, shading="auto", cmap="magma")
plt.plot(ks, exact, "w--", label="exact band")
plt.xlabel("momentum k")
plt.ylabel("energy w")
plt.legend()
plt.show()
