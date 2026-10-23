# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Spectral function A(w) of a half-filled Hubbard chain,
# A(w) = sum_n |<n|c^dagger|GS>|^2 delta(w - E_n + E_0) + hole part.
# It is measured by photoemission and inverse photoemission. For large U
# it splits into two Hubbard bands separated by the Mott gap. We solve
# the N+1 and N-1 electron sectors separately with DMRG.
import numpy as np
from dmrgpy import fermionchain

n = 6 # number of sites
es = np.linspace(-9, 9, 400) # energies

import matplotlib.pyplot as plt
for U in [0.0, 4.0, 8.0]:
    fc = fermionchain.Spinful_Fermionic_Chain_Native(n)
    h = 0
    for i in range(n-1): # hopping
        for cd, c in ((fc.Cdagup, fc.Cup), (fc.Cdagdn, fc.Cdn)):
            h = h - cd[i]*c[i+1] - cd[i+1]*c[i]
    for i in range(n): # Hubbard interaction
        h = h + U*fc.Nup[i]*fc.Ndn[i]
    fc.set_hamiltonian(h)
    fc.maxm = 60 # bond dimension
    fc.nsweeps = 10 # number of DMRG sweeps
    fc.set_conserved_sector(Nf=n, Sz=0) # half filling
    x, A = fc.get_spectral_function(i=n//2, spin="up", nex=10, es=es, delta=0.35)
    plt.plot(x, np.real(A), label="U=%s" % U)
plt.xlabel("energy w")
plt.ylabel("A(w)")
plt.legend()
plt.show()
