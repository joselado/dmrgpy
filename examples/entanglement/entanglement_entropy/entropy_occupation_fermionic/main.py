# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../../src')

# Entanglement entropy of the sites (0,1) with the rest, in the ground
# state of four spinless fermions with first-neighbor hopping. The two
# halves (0,1) and (2,3) are coupled by the hopping c, and we sweep c.
# For c=0 the two halves are decoupled and the entropy is zero.
import numpy as np
from dmrgpy import fermionchain

n = 4 # number of sites
cs = np.linspace(0.0, 1.0, 11) # coupling between the two halves
ss = []
for c in cs:
    fc = fermionchain.Fermionic_Chain(n)
    h = fc.Cdag[0]*fc.C[1] + fc.Cdag[2]*fc.C[3]
    h = h + c*fc.Cdag[1]*fc.C[2] # coupling between the halves
    h = h + h.get_dagger() # add the Hermitian conjugate
    fc.set_hamiltonian(h)
    wf = fc.get_gs() # ground state
    ss.append(wf.get_pair_entropy(0, 1)) # entropy of (0,1) with the rest

import matplotlib.pyplot as plt
plt.plot(cs, ss, marker="o")
plt.xlabel("coupling c between (0,1) and (2,3)")
plt.ylabel("entropy of the pair (0,1)")
plt.show()
