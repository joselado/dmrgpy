# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Quench of the S=1/2 Heisenberg chain. We start from the antiferromagnetic
# ground state of H0 = sum_i (-1)^i Sz_i, and evolve it in time with the
# Heisenberg Hamiltonian H1, so that <Sz_0>(t) = <psi(t)|Sz_0|psi(t)>.
# We compare two time evolution methods, TDVP and the MPO Taylor expansion.
import numpy as np
from dmrgpy import spinchain
from dmrgpy import timedependent

n = 6 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
sc.maxm = 40 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps

h0 = sum([(-1)**i*sc.Sz[i] for i in range(n)])
h1 = 0
for i in range(n-1):
    h1 = h1 + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h0)
wf = sc.get_gs() # initial state
sc.set_hamiltonian(h1) # new Hamiltonian

nt = 300 # number of time steps
dt = 1e-2 # time step
sc.tevol_method = "TDVP"
t1, s1 = timedependent.evolve_and_measure(sc, operator=sc.Sz[0], nt=nt, dt=dt, wf=wf)
sc.tevol_method = "MPO"
t2, s2 = timedependent.evolve_and_measure(sc, operator=sc.Sz[0], nt=nt, dt=dt, wf=wf)
print("Max difference =", np.max(np.abs(s1-s2)))

import matplotlib.pyplot as plt
plt.plot(t2, s2.real, linewidth=4, label="MPO Taylor")
plt.plot(t1, s1.real, "--", linewidth=2, label="TDVP")
plt.xlabel("time")
plt.ylabel("<Sz_0>(t)")
plt.legend()
plt.show()
