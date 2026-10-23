# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Total particle number N = sum_i n_i of a spinless fermion chain with
# hopping t=1 and nearest-neighbor repulsion V=1, as a function of the
# chemical potential mu. N goes down in steps, which are the charge plateaus
# that a gate voltage scan would show.
import numpy as np
from dmrgpy import fermionchain

L = 4 # number of sites

def particle_number(mu):
    fc = fermionchain.Fermionic_Chain(L)
    h = 0
    for i in range(L-1): h = h + fc.Cdag[i]*fc.C[i+1] # hopping
    h = h + h.get_dagger()
    for i in range(L-1): h = h + fc.N[i]*fc.N[i+1] # repulsion V=1
    for i in range(L): h = h + mu*fc.N[i] # chemical potential
    fc.set_hamiltonian(h)
    N = sum(fc.N) # total number operator
    return fc.vev(N).real

mus = np.linspace(-5., 5., 30) # chemical potentials
ns = [particle_number(mu) for mu in mus]

import matplotlib.pyplot as plt
plt.plot(mus, ns, "o-")
plt.xlabel("chemical potential mu")
plt.ylabel("number of electrons")
plt.show()
