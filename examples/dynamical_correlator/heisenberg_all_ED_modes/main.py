# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Dynamical correlator S(w) = sum_n <GS|A|n><n|B|GS> delta(w - E_n + E_0)
# of a spin-1/2 cluster with random exchange, with A=Sy_i and B=Sy_j.
# We compute it with three exact diagonalization methods (submodes),
# which must give the same spectrum. This is what neutron scattering measures.
import numpy as np
from dmrgpy import spinchain

np.random.seed(1) # same random exchange every run
n = 4 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
h = 0
for i in range(n):
    for j in range(n): # random exchange between all pairs
        h = h + sc.Sx[i]*sc.Sx[j]*np.random.random()
        h = h + sc.Sy[i]*sc.Sy[j]*np.random.random()
        h = h + sc.Sz[i]*sc.Sz[j]*np.random.random()
h = h + h.get_dagger()
sc.set_hamiltonian(h)

name = (sc.Sy[0], sc.Sy[1]) # the two operators of the correlator
es = np.linspace(-0.5, 6.0, 1000) # energies
delta = 2e-2 # broadening of the peaks
import matplotlib.pyplot as plt
for submode, style in [("ED", "-"), ("INV", "--"), ("CVM", ":")]:
    x, y = sc.get_dynamical_correlator(mode="ED", name=name, submode=submode,
        es=es, delta=delta)
    plt.plot(x, y.real, style, linewidth=3, label=submode)
plt.xlabel("energy w")
plt.ylabel("S(w)")
plt.legend()
plt.show()
