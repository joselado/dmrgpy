# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Jordan-Wigner transformation. The spin XX chain, H = sum_i Sx_i Sx_{i+1} + Sy_i Sy_{i+1},
# is the same model as free spinless fermions with hopping 1/2.
# The spectral function of Sx_0 in the spin chain must then be equal to
# the one of the fermion operator C_0 (up to a factor 2 in the normalization).
import numpy as np
from dmrgpy import spinchain, fermionchain

n = 6 # number of sites
es = np.linspace(-0.5, 3, 500) # energies
delta = 5e-2 # broadening of the peaks

# spin XX chain
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
sc.set_hamiltonian(h)
(xs, ys) = sc.get_dynamical_correlator(name=(sc.Sx[0], sc.Sx[0]), es=es, delta=delta)

# spinless fermions with hopping 1/2
fc = fermionchain.Fermionic_Chain(n)
h = 0
for i in range(n-1):
    h = h + fc.Cdag[i]*fc.C[i+1]/2.
h = h + h.get_dagger()
fc.set_hamiltonian(h)
(xf, yf) = fc.get_dynamical_correlator(name=(fc.Cdag[0], fc.C[0]), es=es, delta=delta)

print("Max difference =", np.max(np.abs(ys.real-yf.real/2)))

import matplotlib.pyplot as plt
plt.plot(xs, ys.real, linewidth=4, label="spin Sx")
plt.plot(xf, yf.real/2, "--", linewidth=2, label="fermion C")
plt.xlabel("energy w")
plt.ylabel("S(w), first site")
plt.legend()
plt.show()
