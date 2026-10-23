# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Local spin spectrum S(w) of a short Heisenberg chain, computed with the
# correction vector method (CVM) at each frequency. CVM solves the
# resolvent (w - H + i delta)^-1 acting on Sz_0|GS>, so for an exact
# ground state it must agree point by point with exact diagonalization (ED).
import numpy as np
from dmrgpy import spinchain

n = 4 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])

# Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 40 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps

es = np.linspace(-0.5, 4.0, 200) # energies
delta = 0.15 # broadening of the peaks
name = (sc.Sz[0], sc.Sz[0]) # the two operators of the correlator
(e1, s1) = sc.get_dynamical_correlator(name=name, es=es, delta=delta, mode="DMRG", submode="CVM")
(e2, s2) = sc.get_dynamical_correlator(name=name, es=es, delta=delta, mode="ED", submode="ED")
diff = np.max(np.abs(np.array(s1)-np.array(s2)))
print("Max difference CVM-ED =", diff)
assert diff < 1e-6

import matplotlib.pyplot as plt
plt.plot(e2, np.real(s2), linewidth=4, label="ED")
plt.plot(e1, np.real(s1), "--", linewidth=2, label="CVM")
plt.xlabel("energy w")
plt.ylabel("S(w), first site")
plt.legend()
plt.show()
