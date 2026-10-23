# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Local spin spectrum S(w) of the S=1/2 Heisenberg chain with the tangent
# space Krylov method (TaSK), compared with the exact Lehmann sum (ED).
# With maxm large enough to hold the whole Hilbert space TaSK is exact.
# We use the pure Python backend.
import numpy as np
from dmrgpy import spinchain

n = 8 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)], itensor_version="python")

# Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 16 # bond dimension, the whole Hilbert space
sc.nsweeps = 20 # number of DMRG sweeps
sc.cutoff = 0. # keep every Schmidt state

es = np.linspace(0., 3.5, 400) # energies
delta = 0.05 # broadening of the peaks
name = (sc.Sz[n//2], sc.Sz[n//2]) # the two operators of the correlator
(e1, s1) = sc.get_dynamical_correlator(name=name, es=es, delta=delta, submode="TASK")
(e2, s2) = sc.get_dynamical_correlator(name=name, es=es, delta=delta, mode="ED", submode="ED")
diff = np.max(np.abs(np.array(s1)-np.array(s2)))
print("Max difference TaSK-ED =", diff)
assert diff < 1e-6

import matplotlib.pyplot as plt
plt.plot(e2, np.real(s2), linewidth=4, label="ED")
plt.plot(e1, np.real(s1), "--", linewidth=2, label="TaSK")
plt.xlabel("energy w")
plt.ylabel("S(w), middle site")
plt.legend()
plt.show()
