# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Transverse field Ising chain with first and second neighbor couplings,
# H = -G sum sz - J1 sum sx_i sx_{i+1} - J2 sum sx_i sx_{i+2}, in the
# thermodynamic limit (infinite DMRG). The second neighbor term is written
# on a one-site unit cell with get_operator(name, i, group=2), which is the
# cell two steps to the right. The same model on a two-site cell must agree.
import numpy as np
from dmrgpy import infinitechain

G, J1 = 2.5, 1.0 # transverse field and first neighbor coupling

def energy(n_uc, J2):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"]*n_uc)
    ic.maxm, ic.maxiter, ic.etol = 8, 300, 1e-10 # bond dimension and accuracy
    h = 0
    for i in range(n_uc):
        h = h - 2*G*ic.SzC[i] # sigma^z = 2 Sz
        for r, J in ((1, J1), (2, J2)): # neighbor at distance r
            k = i + r
            if k < n_uc: other = ic.SxC[k]
            else: other = ic.get_operator("Sx", k % n_uc, group=k // n_uc)
            h = h - 4*J*ic.SxC[i]*other # sigma^x = 2 Sx
    ic.set_hamiltonian(h)
    return ic.gs_energy()

j2s = np.linspace(0., 1., 6)
e1 = np.array([energy(1, J2) for J2 in j2s]) # one-site cell, reach 2
e2 = np.array([energy(2, J2) for J2 in j2s]) # two-site cell, reach 1
print("Max difference =", np.max(np.abs(e1-e2)))

import matplotlib.pyplot as plt
plt.plot(j2s, e1, linewidth=4, label="one-site cell")
plt.plot(j2s, e2, "--", linewidth=2, label="two-site cell")
plt.xlabel("J2")
plt.ylabel("energy per site")
plt.legend()
plt.show()
