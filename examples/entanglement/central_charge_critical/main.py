# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Central charge of the critical transverse field Ising chain. The bond
# entropy of a critical open chain follows S(l) = c/6 log(2L/pi sin(pi l/L)) + const,
# so the fit to the entropy profile gives the central charge c (c=1/2 for Ising).
import numpy as np
from dmrgpy import spinchain

n = 100 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
h = 0
for i in range(n-1): h = h + sc.Sz[i]*sc.Sz[i+1] # Ising coupling
for i in range(n): h = h + 0.5*sc.Sx[i] # transverse field
sc.set_hamiltonian(h)
sc.maxm = 200 # bond dimension, large for a critical system
wf = sc.get_gs()
c = wf.get_CFT_central_charge() # central charge from the fit
print("Central charge", c)

ls = np.arange(1, n) # bonds
s = [wf.get_bond_entropy(l-1, l) for l in ls] # entropy of each bond
fit = c/6*np.log(2*n/np.pi*np.sin(np.pi*ls/n))
fit = fit + np.mean(s - fit) # constant of the fit

import matplotlib.pyplot as plt
plt.plot(ls, s, "o", label="DMRG")
plt.plot(ls, fit, "r-", label="CFT fit, c=%.3f" % c)
plt.xlabel("bond")
plt.ylabel("Entanglement entropy")
plt.legend()
plt.show()
