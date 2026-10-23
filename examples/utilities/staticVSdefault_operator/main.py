# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Repeated application of the Hamiltonian, |psi> -> H|psi>, on a random
# MPS of a Heisenberg chain. The norm <psi|psi> grows as H^(2k), and it
# is a check that a StaticOperator, built once, gives the same MPS as
# the default operator, which is rebuilt on every call.
import numpy as np
from dmrgpy import spinchain
from dmrgpy.multioperatortk.staticoperator import StaticOperator

n = 20 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
wf0 = sc.random_mps()
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sh = StaticOperator(h, sc) # static version of h

def norms(O):
    wf = wf0.copy() ; ds = []
    for i in range(10):
        ds.append(wf.dot(wf).real)
        wf = O*wf
    return np.array(ds)

d1 = norms(h)
d2 = norms(sh)
print("Max relative difference =", np.max(np.abs(d1-d2)/d1))

import matplotlib.pyplot as plt
plt.semilogy(d1, linewidth=4, label="default")
plt.semilogy(d2, "--", linewidth=2, label="static")
plt.xlabel("iteration k")
plt.ylabel("<psi|psi>")
plt.legend()
plt.show()
