# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Spectral function of a Kondo impurity, a Hubbard site coupled to a
# logarithmic bath of 2*nb+1 levels (a Wilson-like discretization).
# A(w) = sum_n |<n|c_up|GS>|^2 delta(w - E_n + E_0) is what an STM
# measures on top of the impurity. We compute it with DMRG and KPM.
import numpy as np
from dmrgpy import fermionchain

nb = 6 # bath levels on each side
n = 2*nb + 2 # total number of sites
ic = nb # impurity site
fc = fermionchain.Spinful_Fermionic_Chain(n)
U = 4.0 # Hubbard interaction
h = U*(fc.Nup[ic]-0.5)*(fc.Ndn[ic]-0.5)
eps = np.array([0.4**(nb-i-1) for i in range(nb)]) # log-spaced levels
eps = np.sort(np.concatenate([-eps, [0.], eps])) # bath energies
inds = list(range(nb)) + [nb+1] + [nb+i+2 for i in range(nb)] # bath sites
t = 1./nb # coupling to the impurity
for e, i in zip(eps, inds):
    h = h + t*fc.Cdagup[i]*fc.Cup[ic] + t*fc.Cdagdn[i]*fc.Cdn[ic]
    h = h + e*(fc.Nup[i] + fc.Ndn[i])
h = h + h.get_dagger()
fc.set_hamiltonian(h)
fc.maxm = 20 # bond dimension
fc.kpmmaxm = 20 # bond dimension in KPM

name = (fc.Cdagup[ic], fc.Cup[ic])
es = np.linspace(-3.0, 3.0, 300) # energies
x, y = fc.get_dynamical_correlator(name=name, es=es, delta=0.1)

import matplotlib.pyplot as plt
plt.plot(x, y.real)
plt.xlabel("energy w")
plt.ylabel("A(w), impurity")
plt.show()
