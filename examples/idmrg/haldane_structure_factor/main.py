# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Magnon dispersion of the Haldane phase, the S=1 Heisenberg chain
# H = sum_j S_j.S_{j+1}, in the thermodynamic limit (iDMRG/VUMPS).
# spectral_weights gives the energies and the weights of the excitations
# that Sz_k creates. The lowest branch has a gap of 0.41 at k=pi,
# which inelastic neutron scattering measures.
import numpy as np
from dmrgpy import infinitechain

ic = infinitechain.Infinite_Spin_Chain(["1"], itensor_version="python")
Sc = [ic.SxC[0], ic.SyC[0], ic.SzC[0]] # operators on the unit cell
Sr = [ic.SxR[0], ic.SyR[0], ic.SzR[0]] # operators on the next cell
ic.set_hamiltonian(Sc[0]*Sr[0] + Sc[1]*Sr[1] + Sc[2]*Sr[2])
ic.gs_method = "vumps"
ic.maxm = 10 # bond dimension
print("Energy per site =", ic.gs_energy(), "(exact -1.401484)")

ks = np.linspace(0.0, 2*np.pi, 25)
band, weight = [], []
for k in ks:
    e, w, tot = ic.spectral_weights("Sz", k, n=12, return_total=True)
    band.append(e[0]) # lowest branch
    weight.append(tot) # total weight of Sz_k
print("Gap at k=pi =", band[len(ks)//2], "(exact 0.4105)")

import matplotlib.pyplot as plt
plt.plot(ks, band, "o-", label="lowest branch")
plt.plot(ks, weight, "s--", label="S(k), total weight")
plt.xlabel("momentum k")
plt.ylabel("energy / weight")
plt.legend()
plt.show()
