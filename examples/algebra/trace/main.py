# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Trace of an operator, Tr(O) = sum_s <s|O|s>, summed over all the states
# of the chain, and the trace of its inverse. Here O = H^2 + i, with H the
# Heisenberg Hamiltonian. The inverse needs no diagonalization of H.
# The trace gives, for example, the infinite temperature average.
# We compute it with matrix product states (MPS) and compare with ED.
import numpy as np
from dmrgpy import spinchain

ns = [2, 4, 6] # chain lengths, ED cost grows exponentially
tr_ed, tr_mps = [], []
for n in ns:
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
    # Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    o = h*h + 1j
    tr_ed.append(sc.trace(o, mode="ED"))
    tr_mps.append(sc.trace(o, mode="DMRG"))
    itr_ed = sc.inverse_trace(o, mode="ED")
    itr_mps = sc.inverse_trace(o, mode="DMRG")
    print("n =", n, " inverse trace difference =", abs(itr_ed-itr_mps))
print("Max difference =", np.max(np.abs(np.array(tr_ed)-tr_mps)))

import matplotlib.pyplot as plt
plt.plot(ns, np.abs(tr_ed), linewidth=4, label="ED")
plt.plot(ns, np.abs(tr_mps), "--", linewidth=2, label="MPS")
plt.yscale("log")
plt.xlabel("number of sites n")
plt.ylabel("|Tr(H^2 + i)|")
plt.legend()
plt.show()
