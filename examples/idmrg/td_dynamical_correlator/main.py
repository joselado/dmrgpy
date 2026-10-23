# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Dynamical spin correlator S(k,w) of an infinite chain in the paramagnetic
# phase of a transverse-field Ising model, H = 1.4 Sz + Sx_i Sx_{i+1}.
# The ground state comes from iDMRG. The dynamics comes from real-time
# TDVP on a finite window capped with the iDMRG environments.
# The bright line is the magnon band, seen in neutron scattering.
import numpy as np
from dmrgpy import infinitechain

ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
ic.gs_method = "idmrg"
ic.set_hamiltonian(1.4*ic.SzC[0] + ic.SxC[0]*ic.SxR[0])
ic.maxm, ic.maxiter, ic.etol, ic.niter = 20, 300, 1e-12, 150 # accuracy
print("Energy density =", ic.gs_energy())

ks = np.linspace(-np.pi, np.pi, 31) # momenta
es = np.linspace(-1, 6, 100) # energies
(_, _, Skw) = ic.td_dynamical_correlator("Sx", 0, "Sx", n_window=14, dt=0.1,
        nt=60, maxdim=60, cutoff=1e-10, niter=50, x_values=range(-4, 5),
        ks=ks, es=es, delta=0.15)

import matplotlib.pyplot as plt
plt.pcolormesh(ks, es, np.abs(Skw).T, shading="auto", cmap="magma")
plt.colorbar(label="|S(k,w)|")
plt.xlabel("momentum k")
plt.ylabel("energy w")
plt.show()
