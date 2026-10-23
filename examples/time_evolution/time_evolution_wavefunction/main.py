# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Time evolution of a wavefunction, |psi(t)> = exp(-i H t)|psi(0)>.
# We start from the Neel state of four spins, which is the ground state of a
# staggered field, and evolve it with the Heisenberg coupling.
# The entanglement entropy of the first site measures how the two spins entangle.
import numpy as np
from dmrgpy import spinchain
from dmrgpy.timeevolution import imaginary_exponential as evolve

n = 4 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])

h0 = 0 # staggered field, to create the initial state
for i in range(n):
    h0 = h0 + (-1)**i*sc.Sz[i]
sc.set_hamiltonian(h0)
wf0 = sc.get_gs() # initial state

h = 0 # Heisenberg Hamiltonian, to evolve the state
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
h = h + h.get_dagger()

ts = np.linspace(0., 4., 40) # times
ss = [evolve(h, wf0, ts=[t])[0].get_site_entropy(0) for t in ts]

import matplotlib.pyplot as plt
plt.plot(ts, ss, "o-")
plt.xlabel("time")
plt.ylabel("entropy of the first site")
plt.show()
