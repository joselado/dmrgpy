# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Many-body states and operators behave as vectors and matrices. We
# build a random state |2>, multiply it by a complex number, sum states,
# apply the Heisenberg Hamiltonian H and compute overlaps <2|H|2>.
# The expectation value grows linearly with the number of sites.
import numpy as np
from dmrgpy import spinchain

def overlap(n):
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    psi1 = sc.random_state(mode="MPS") # random state |1>
    psi2 = sc.random_state(mode="MPS") # random state |2>
    psi3 = (0.3+0.1j)*psi2 # |3> = z |2>
    psi4 = h*psi2 # |4> = H |2>
    psi5 = psi2 + psi4 # |5> = |2> + |4>
    return psi2.dot(psi4) # overlap <2|H|2>

print("Overlap <2|H|2> for n = 6:", overlap(6))
ns = range(4, 13, 2) # chain lengths
out = [overlap(n).real for n in ns]

import matplotlib.pyplot as plt
plt.plot(ns, out, "o-")
plt.xlabel("number of sites n")
plt.ylabel("<2|H|2>")
plt.show()
