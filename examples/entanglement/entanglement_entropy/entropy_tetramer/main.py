# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../../src')

# Entanglement entropy of the pair of sites (0,1) with the pair (2,3) in a
# Heisenberg chain of four sites, S = -Tr rho log(rho), where rho is the
# reduced density matrix of the pair. We change the coupling in the middle.
# At zero coupling the two dimers are decoupled and the entropy vanishes.
import numpy as np
from dmrgpy import spinchain

n = 4 # number of sites
cs = np.linspace(0, 1, 20) # coupling between the two pairs
ss = []
for c01 in cs:
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
    # Heisenberg Hamiltonian, the middle bond has coupling c01
    h = 0
    for i in range(n-1):
        c = 1.0
        if i == 1: c = c01
        h = h + c*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
    sc.set_hamiltonian(h)
    wf = sc.get_gs() # ground state
    ss.append(wf.get_pair_entropy(0, 1)) # entropy of sites (0,1) with the rest

import matplotlib.pyplot as plt
plt.plot(cs, ss, marker="o")
plt.xlabel("coupling between the two pairs")
plt.ylabel("entanglement entropy of the pair (0,1)")
plt.show()
