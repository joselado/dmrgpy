# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Quench of a spin chain. We take the ground state of a staggered field,
# switch it off and let the chain evolve with the Heisenberg Hamiltonian.
# We follow <Sz_0>(t) with TDVP and compare it with exact diagonalization.
import numpy as np
from dmrgpy import spinchain
from dmrgpy import timedependent

n = 4 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
h0, h1 = 0, 0
for i in range(n):
    h0 = h0 + (-1)**i*sc.Sz[i] # staggered field
for i in range(n-1):
    h1 = h1 + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h0)
sc.maxm = 40 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps
wf = sc.get_gs() # DMRG ground state
wfED = sc.get_gs(mode="ED") # ED ground state
sc.set_hamiltonian(h1) # quench to the Heisenberg chain

nt, dt = 50, 0.05 # number of time steps and time step
(ts, sz) = timedependent.evolve_and_measure(sc, operator=sc.Sz[0], nt=nt, dt=dt, wf=wf)
(ts, szED) = timedependent.evolve_and_measure(sc, operator=sc.Sz[0], nt=nt, dt=dt, wf=wfED, mode="ED")
diff = np.max(np.abs(sz-szED))
print("Max difference TDVP-ED =", diff)
assert diff < 1e-4

import matplotlib.pyplot as plt
plt.plot(ts, szED.real, linewidth=4, label="ED")
plt.plot(ts, sz.real, "--", linewidth=2, label="TDVP")
plt.xlabel("time")
plt.ylabel("<Sz_0>")
plt.legend()
plt.show()
