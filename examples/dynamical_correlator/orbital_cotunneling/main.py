# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Orbital cotunneling in a two-orbital system. The electron hops from
# orbital 1 to orbital 0 (and back) with the operator A = Cdag_0 C_1, and
# the spectrum S(w) = sum_n |<n|A|GS>|^2 delta(w - E_n + E_0) gives the
# energy of the excitations that this process creates.
import numpy as np
from dmrgpy import fermionchain

n = 2 # number of orbitals
fc = fermionchain.Spinful_Fermionic_Chain(n)

tij = [[0.5, 1.0], [1.0, -0.5]] # hopping matrix
h = 0
for i in range(n):
    for j in range(n):
        h = h + tij[i][j]*fc.Cdagup[i]*fc.Cup[j]
        h = h + tij[i][j]*fc.Cdagdn[i]*fc.Cdn[j]
fc.set_hamiltonian(h)

wf = fc.get_gs(mode="ED") # ground state
A01 = fc.Cdagup[0]*fc.Cup[1] # hop from 1 to 0
A10 = fc.Cdagup[1]*fc.Cup[0] # hop from 0 to 1
A01 = A01 - wf.dot(A01*wf) # remove the expectation value
A10 = A10 - wf.dot(A10*wf)

es = np.linspace(0, 4, 500) # energies
es1, d01 = fc.get_dynamical_correlator(name=(A01, A01.get_dagger()), mode="ED", submode="INV", es=es)
es2, d10 = fc.get_dynamical_correlator(name=(A10, A10.get_dagger()), mode="ED", submode="INV", es=es)

import matplotlib.pyplot as plt
plt.plot(es1, d01.real, label="1 -> 0")
plt.plot(es2, d10.real, label="0 -> 1")
plt.xlabel("energy w")
plt.ylabel("S(w)")
plt.legend()
plt.show()
