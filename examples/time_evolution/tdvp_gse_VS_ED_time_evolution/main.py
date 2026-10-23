# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Quench dynamics of <Sz_0>(t). We start from the ground state of a
# staggered field, and evolve with the Heisenberg Hamiltonian. The evolution
# uses one-site TDVP with global subspace expansion (TDVP_GSE), and we
# compare with exact diagonalization (ED).
import numpy as np
from dmrgpy import spinchain, timedependent

n = 4 # number of sites
sc = spinchain.Spin_Chain([2 for i in range(n)]) # spin 1/2
sc.tevol_method = "TDVP_GSE"
h0 = 0 # initial Hamiltonian
h1 = 0 # Hamiltonian for the evolution
for i in range(n):
    h0 = h0 + (-1)**i*sc.Sz[i] # staggered field
for i in range(n-1):
    h0 = h0 + 0.3*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1])
    h1 = h1 + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h0)
wf = sc.get_gs() # DMRG ground state
wf_ed = sc.get_gs(mode="ED") # ED ground state
sc.set_hamiltonian(h1) # quench

(t, sz) = timedependent.evolve_and_measure(sc, operator=sc.Sz[0], nt=50, dt=0.05, wf=wf)
(t, sz_ed) = timedependent.evolve_and_measure(sc, operator=sc.Sz[0], nt=50, dt=0.05, wf=wf_ed, mode="ED")
diff = np.max(np.abs(sz - sz_ed))
print("Max difference =", diff)
assert diff < 1e-4

import matplotlib.pyplot as plt
plt.plot(t, sz_ed.real, linewidth=4, label="ED")
plt.plot(t, sz.real, "--", linewidth=2, label="TDVP_GSE")
plt.xlabel("time")
plt.ylabel("<Sz_0>(t)")
plt.legend()
plt.show()
