# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Charge response of a free fermionic chain, H = sum_ij m_ij c^dagger_i c_j.
# The charge correlator <dN_i(t) dN_j> has a Lehmann spectrum made of
# particle-hole excitations between occupied and empty levels (Lindhard
# form). We compute it with DMRG and compare with that formula.
import numpy as np
import scipy.linalg as lg
from dmrgpy import fermionchain

n = 6 # number of sites
fc = fermionchain.Fermionic_Chain(n)
np.random.seed(1)
m = np.random.random((n,n)) + 1j*np.random.random((n,n))
m = m + np.conjugate(m).T # hopping matrix, Hermitian
h = 0
for i in range(n):
    for j in range(n):
        h = h + m[i,j]*fc.Cdag[i]*fc.C[j]
fc.set_hamiltonian(h)
fc.maxm = 40 # bond dimension
fc.nsweeps = 10 # number of DMRG sweeps

i = 1 # site of the density operators
ni = fc.vev(fc.N[i])
name = (fc.N[i]-ni*fc.Id, fc.N[i]-ni*fc.Id) # density fluctuation
es = np.linspace(-0.5, 6.0, 500) # energies
delta = 0.1 # broadening
x, y = fc.get_dynamical_correlator(mode="DMRG", name=name, es=es, delta=delta)

(ei, vi) = lg.eigh(m)
eta = 0.06 # broadening of the formula
y2 = np.zeros(len(es)) # Lindhard formula
for a in range(n):
    for b in range(n):
        if ei[a] < 0.0 and ei[b] > 0.0: # occupied a, empty b
            w = abs(vi[i,a])**2*abs(vi[i,b])**2
            y2 += w*eta/(eta**2 + (es+ei[a]-ei[b])**2)/np.pi

import matplotlib.pyplot as plt
plt.plot(x, y2, linewidth=4, label="Lindhard")
plt.plot(x, y.real, "--", linewidth=2, label="DMRG")
plt.xlabel("energy w")
plt.ylabel("charge correlator")
plt.legend()
plt.show()
