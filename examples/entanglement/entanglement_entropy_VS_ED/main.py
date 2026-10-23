# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Entanglement of the ground state of a 4-site Heisenberg chain with a
# small field on the first site. We compute the entropy of one site, the
# entropy of a pair of sites and the mutual information between two
# sites, S = -Tr(rho log rho) with rho the reduced density matrix.
# We compare DMRG with exact diagonalization (ED).
import numpy as np
from dmrgpy import spinchain

n = 4 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
h = 0.3*sc.Sz[0] # small field
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 20 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps

def entropies(wf):
    return [sc.get_site_entropy(wf, 0), sc.get_pair_entropy(wf, 0, 1),
            sc.get_mutual_information(wf, 0, 2)]
s_ed = entropies(sc.get_gs(mode="ED"))
s_dmrg = entropies(sc.get_gs(mode="DMRG"))
print("Max difference DMRG-ED =", np.max(np.abs(np.array(s_dmrg) - s_ed)))

import matplotlib.pyplot as plt
x = np.arange(3)
plt.bar(x-0.2, s_ed, 0.4, label="ED")
plt.bar(x+0.2, s_dmrg, 0.4, label="DMRG")
plt.xticks(x, ["site 0", "pair 0,1", "mutual 0,2"])
plt.ylabel("entropy")
plt.legend()
plt.show()
