# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Local spin spectrum S(w) = sum_n |<n|Sx_0|GS>|^2 delta(w - E_n + E_0) of
# the S=1/2 Heisenberg chain, computed with the default C++ DMRG backend
# and with the Julia backend (ITensors.jl). Both must give the same peaks,
# which are the spin excitations seen in neutron scattering.
import numpy as np
from dmrgpy import spinchain

def get_chain(n):
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
    # Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = 40 # bond dimension
    sc.nsweeps = 10 # number of DMRG sweeps
    return sc

n = 8 # number of sites
es = np.linspace(-0.5, 3.0, 300) # energies
sc = get_chain(n) # C++ backend
(e1, s1) = sc.get_dynamical_correlator(name=(sc.Sx[0], sc.Sx[0]), es=es, delta=0.05)
sc = get_chain(n)
sc.setup_julia() # Julia backend
(e2, s2) = sc.get_dynamical_correlator(name=(sc.Sx[0], sc.Sx[0]), es=es, delta=0.05)
print("Max difference =", np.max(np.abs(s1-s2)))

import matplotlib.pyplot as plt
plt.plot(e1, s1.real, linewidth=4, label="C++")
plt.plot(e2, s2.real, "--", linewidth=2, label="Julia")
plt.xlabel("energy w")
plt.ylabel("S(w), first site")
plt.legend()
plt.show()
