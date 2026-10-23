# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../../src')

# Entanglement entropy of the first site of a Heisenberg dimer,
# S = -Tr(rho log rho), with rho the reduced density matrix of the site.
# The singlet ground state has the maximum entropy, log(2). A local
# field Bz on the first site polarizes it and the entropy goes to zero.
import numpy as np
from dmrgpy import spinchain

sc = spinchain.Spin_Chain(["S=1/2" for i in range(2)])
h = sc.Sx[0]*sc.Sx[1] + sc.Sy[0]*sc.Sy[1] + sc.Sz[0]*sc.Sz[1]

bzs = np.linspace(0., 4.0, 20) # local field on the first site
ss = [] # entropies
for bz in bzs:
    sc.set_hamiltonian(h + bz*sc.Sz[0])
    wf = sc.get_gs() # ground state
    ss.append(sc.get_site_entropy(wf, 0)) # entropy of the first site

import matplotlib.pyplot as plt
plt.plot(bzs, ss, "o-")
plt.axhline(np.log(2), color="gray", linestyle="--")
plt.xlabel("local field Bz on site 0")
plt.ylabel("entanglement entropy")
plt.show()
