# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Non-Hermitian chain of spinless fermions with gain and loss, in a ring of
# 8 sites whose last bond has strength c. For c=0 the chain is open and
# has edge zero modes, so the ground state is degenerate. Closing the
# ring with c splits the levels. We use exact diagonalization (ED).
import numpy as np
from dmrgpy import fermionchain

def get_spectrum(c, eta=1.0):
    ns = 8 # number of sites
    t = np.zeros((ns, ns), dtype=complex) # hopping matrix
    for i in range(ns-1):
        t[i, i+1] = t[i+1, i] = 1.0
    t[0, ns-1] = t[ns-1, 0] = c # closing bond
    for i in range(ns): # gain and loss pattern
        t[i, i] = 1j*eta*[1, -1, -1, 1][i % 4]
    fc = fermionchain.Fermionic_Chain(ns)
    h = 0
    for i in range(ns):
        for j in range(ns):
            h = h + t[i, j]*fc.Cdag[i]*fc.C[j]
    fc.set_hamiltonian(h)
    (es, wfs) = fc.get_excited_states(mode="ED", n=10)
    return es - es[0], fc.get_gs_degeneracy(mode="ED", delta=0.1)

cs = np.linspace(0.0, 1.0, 20) # strength of the closing bond
xs, ys = [], []
for c in cs:
    (es, deg) = get_spectrum(c)
    print("c =", round(c, 2), " degeneracy =", deg)

    xs += [c]*len(es)
    ys += list(es.real)

import matplotlib.pyplot as plt
plt.scatter(xs, ys, s=10)
plt.xlabel("closing bond c")
plt.ylabel("Re(E-E0)")
plt.show()
