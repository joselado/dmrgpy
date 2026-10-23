# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Quench dynamics <Sz_0(t)> of a Heisenberg chain. We start from the ground
# state of a staggered field (Neel order) and let it evolve with the
# Heisenberg Hamiltonian. We use time evolving block decimation (TEBD),
# available in the pure Python backend, and compare with exact diagonalization.
import numpy as np
from dmrgpy import spinchain, timedependent

n = 4 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
sc.setup_python() # TEBD is in the pure Python backend
sc.tevol_method = "TEBD"

h0 = sum((-1)**i*sc.Sz[i] for i in range(n)) # staggered field
h1 = 0
for i in range(n-1): # Heisenberg Hamiltonian
    h1 = h1 + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h0)
wf = sc.get_gs() # initial state
wfed = sc.get_gs(mode="ED")
sc.set_hamiltonian(h1) # quench

nt, dt = 50, 0.05 # number of time steps and time step
(t1, s1) = timedependent.evolve_and_measure(sc, operator=sc.Sz[0], nt=nt, dt=dt, wf=wf)
(t2, s2) = timedependent.evolve_and_measure(sc, operator=sc.Sz[0], nt=nt, dt=dt, wf=wfed, mode="ED")
diff = np.max(np.abs(s1-s2))
print("Max difference TEBD-ED =", diff)
assert diff < 1e-4

import matplotlib.pyplot as plt
plt.plot(t2, s2.real, linewidth=4, label="ED")
plt.plot(t1, s1.real, "--", linewidth=2, label="TEBD")
plt.xlabel("time")
plt.ylabel("<Sz_0>(t)")
plt.legend()
plt.show()
