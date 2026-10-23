# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Quench of a spin-1/2 Heisenberg chain starting from a product state,
# the ground state of a staggered field. We follow <Sz_0>(t) with one-site
# TDVP plus global subspace expansion (GSE), which grows the bond dimension,
# and compare with exact diagonalization (ED). Uses the Julia backend.
import numpy as np
from dmrgpy import spinchain
from dmrgpy import timedependent

n = 6 # number of sites
nt = 40 # number of time steps
dt = 0.05 # time step
sc = spinchain.Spin_Chain([2 for i in range(n)])
sc.setup_julia()
sc.tevol_method = "TDVP_GSE"
sc.tdvp_gse_sweeps = nt # steps with subspace expansion

h0 = 0 # staggered field, product ground state
for i in range(n): h0 = h0 + (-1)**i*sc.Sz[i]
h1 = 0 # Heisenberg chain after the quench
for i in range(n-1):
    h1 = h1 + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h0)
wf = sc.get_gs() # DMRG state
wfed = sc.get_gs(mode="ED") # ED state
sc.set_hamiltonian(h1) # quench

op = sc.Sz[0]
t, sz = timedependent.evolve_and_measure(sc, operator=op, nt=nt, dt=dt, wf=wf)
t, szed = timedependent.evolve_and_measure(sc, operator=op, nt=nt, dt=dt,
    wf=wfed, mode="ED")
diff = np.max(np.abs(sz - szed))
print("Max difference =", diff)
assert diff < 1e-4

import matplotlib.pyplot as plt
plt.plot(t, szed.real, linewidth=4, label="ED")
plt.plot(t, sz.real, "--", linewidth=2, label="TDVP_GSE")
plt.xlabel("time")
plt.ylabel("<Sz_0>")
plt.legend()
plt.show()
