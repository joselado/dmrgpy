# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Correlator <C^dagger_0 C_j> in a spinful fermionic chain. The operator
# has a long Jordan-Wigner string, so for big systems we first build it as
# a product of MPOs with toMPO. We compare it with the direct operator.
import numpy as np
from dmrgpy import fermionchain
from dmrgpy.multioperatortk.longoperator import toMPO

ns = 10 # number of spinful sites
fc = fermionchain.Spinful_Fermionic_Chain(ns)

# hopping t=1 and chemical potential mu
mu = -1.0
h = 0
for i in range(ns-1):
    h = h + fc.Cdagup[i]*fc.Cup[i+1] + fc.Cdagdn[i]*fc.Cdn[i+1]
h = h + h.get_dagger()
for i in range(ns):
    h = h + mu*(fc.Cdagup[i]*fc.Cup[i] + fc.Cdagdn[i]*fc.Cdn[i])
fc.set_hamiltonian(h)
fc.maxm = 40 # bond dimension
fc.nsweeps = 10 # number of DMRG sweeps
wf = fc.get_gs() # ground state

full, prod = [], []
for j in range(ns):
    op = fc.Cdagup[0]*fc.Cup[j]
    full.append(wf.dot(op*wf).real) # direct operator
    prod.append(wf.dot(toMPO(fc,op)*wf).real) # product of MPOs
print("Max difference =", np.max(np.abs(np.array(full)-prod)))

import matplotlib.pyplot as plt
plt.plot(range(ns), full, linewidth=4, label="direct")
plt.plot(range(ns), prod, "--", linewidth=2, label="product of MPOs")
plt.xlabel("site j")
plt.ylabel("<C^dagger_0 C_j>")
plt.legend()
plt.show()
