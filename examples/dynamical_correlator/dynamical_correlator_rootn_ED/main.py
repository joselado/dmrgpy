# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Spin spectrum S(w) of the Heisenberg chain from a correction vector.
# The root-N method splits one big Krylov space of dimension nkry into N
# smaller steps. At a fixed small nkry, larger N should be closer to the
# exact Lehmann sum (ED), mainly in the upper part of the spectrum.
import numpy as np
from dmrgpy import spinchain

n = 6 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])

# Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)

name = (sc.Sz[0], sc.Sz[0]) # the two operators of the correlator
es = np.linspace(2.0, 4.5, 40) # upper part of the spectrum
delta = 0.2 # broadening of the peaks
nkry = 4 # small Krylov dimension, to see the error

def spectrum(submode, **kwargs):
    return sc.get_dynamical_correlator(mode="ED", submode=submode, name=name,
            es=es, delta=delta, **kwargs)

(x, yexact) = spectrum("ED")
import matplotlib.pyplot as plt
plt.plot(x, yexact.real, "k", linewidth=4, label="exact")
for N in [1, 4, 8]: # N=1 is the usual correction vector
    (x, y) = spectrum("ROOTN", N=N, nkry=nkry)
    print("N =", N, " max error =", np.max(np.abs(y - yexact)))
    plt.plot(x, y.real, "--", linewidth=2, label="root-N, N=%d" % N)
plt.xlabel("energy w")
plt.ylabel("S(w), first site")
plt.legend()
plt.show()
