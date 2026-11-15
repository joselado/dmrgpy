# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Lowest energy levels of the S=1/2 Heisenberg chain, from the Krylov
# space of two imaginary-time evolved random states (arXiv:2509.07709),
# on the Julia backend. We diagonalize G = -H in the space of the
# snapshots, so the energies are E = -lambda, lowest first. The gaps
# E_k - E_0 are what inelastic neutron scattering measures.
# We check the levels against exact diagonalization (ED).
import numpy as np
from dmrgpy import spinchain

n = 6 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
sc.setup_julia() # ITensors.jl backend
sc.maxm = 40 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps

# Heisenberg Hamiltonian in a field, H = sum_i S_i . S_{i+1} + B sum_i Sz_i
# Two random states find at most two copies of a degenerate level, so a
# small field splits the triplets and each level appears once.
B = 0.2 # magnetic field, below the spin gap
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
for i in range(n): h = h + B*sc.Sz[i]
sc.set_hamiltonian(h)

nlev = 4 # number of levels we compare
ev = sc.get_krylov_spectrum(nk=20, dt=0.5)[0] # 20 snapshots per state
eved = sc.get_krylov_spectrum(mode="ED")[0] # every eigenvalue, same order
e = (-ev).real[:nlev] # energies, lowest first
eed = (-eved).real[:nlev]
print("Max difference Julia-ED =", np.max(np.abs(e-eed)))
assert np.max(np.abs(e-eed)) < 1e-3

import matplotlib.pyplot as plt
plt.plot(range(nlev), eed, "o", markersize=14, label="ED")
plt.plot(range(nlev), e, "x--", markersize=10, label="Julia, Krylov")
plt.xlabel("level index k")
plt.ylabel("energy E_k")
plt.legend()
plt.show()
