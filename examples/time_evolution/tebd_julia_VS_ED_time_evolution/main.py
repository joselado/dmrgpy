# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Quench of the S=1/2 Heisenberg chain. The chain starts in the Neel state,
# the ground state of a staggered field, and evolves with the Heisenberg
# Hamiltonian, <Sz_0(t)> = <psi(t)|Sz_0|psi(t)>. We use TEBD (Trotter
# gates) on the Julia backend and compare with exact diagonalization (ED).
import numpy as np
from dmrgpy import spinchain, timedependent

n = 4 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
sc.setup_julia() # Julia backend
sc.tevol_method = "TEBD"
h0 = 0 # staggered field
for i in range(n): h0 = h0 + (-1)**i*sc.Sz[i]
h1 = 0 # Heisenberg Hamiltonian
for i in range(n-1):
    h1 = h1 + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h0)
wf, wfED = sc.get_gs(), sc.get_gs(mode="ED") # Neel states
sc.set_hamiltonian(h1)

nt, dt = 50, 0.05 # number of time steps, time step
(ts, sz) = timedependent.evolve_and_measure(sc, operator=sc.Sz[0], nt=nt, dt=dt, wf=wf)
(tsED, szED) = timedependent.evolve_and_measure(sc, operator=sc.Sz[0], nt=nt, dt=dt,
                                                wf=wfED, mode="ED")
print("Max difference =", np.max(np.abs(sz-szED)))
assert np.max(np.abs(sz-szED)) < 1e-4

import matplotlib.pyplot as plt
plt.plot(tsED, szED.real, linewidth=4, label="ED")
plt.plot(ts, sz.real, "--", linewidth=2, label="TEBD (Julia)")
plt.xlabel("time t")
plt.ylabel("<Sz_0>(t)")
plt.legend()
plt.show()
