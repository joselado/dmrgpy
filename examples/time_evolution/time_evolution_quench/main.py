# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Quantum quench in the S=1/2 chain. We start from the ground state of a
# staggered field, H0 = sum_i (-1)^i Sz_i, and evolve it in time with
# the Heisenberg Hamiltonian H = sum_i S_i . S_{i+1}. We follow <Sz_0>(t)
# with DMRG and with exact diagonalization (ED).
import numpy as np
from dmrgpy import spinchain, timedependent

n = 6 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
h0 = 0 # initial Hamiltonian
for i in range(n):
    h0 = h0 + (-1)**i*sc.Sz[i]
h1 = 0 # Heisenberg Hamiltonian
for i in range(n-1):
    h1 = h1 + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.maxm = 40 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps

sc.set_hamiltonian(h0)
wf1 = sc.get_gs() # initial state, DMRG
wf2 = sc.get_gs(mode="ED") # initial state, ED
sc.set_hamiltonian(h1) # quench

kw = dict(operator=sc.Sz[0], nt=1000, dt=0.01)
(t1, s1) = timedependent.evolve_and_measure(sc, wf=wf1, **kw)
(t2, s2) = timedependent.evolve_and_measure(sc, wf=wf2, mode="ED", **kw)
print("Max difference DMRG-ED =", np.max(np.abs(np.array(s1)-np.array(s2))))

import matplotlib.pyplot as plt
plt.plot(t2, np.real(s2), linewidth=4, label="ED")
plt.plot(t1, np.real(s1), "--", linewidth=2, label="DMRG")
plt.xlabel("time t")
plt.ylabel("<Sz_0>(t)")
plt.legend()
plt.show()
