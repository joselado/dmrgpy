# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Products of operators as MPOs. We convert a non-Hermitian operator A
# (Heisenberg chain plus imaginary terms) with toMPO, multiply it with its
# adjoint, and compute the traces. The traces grow with the Hilbert space.
from dmrgpy import spinchain

ns = [2,4,6,8] # chain lengths
tr = {"A":[], "A*A":[], "A*A^dagger":[], "A*A^dagger*A":[]}
for n in ns:
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
    A0 = 0.2j # constant imaginary term
    for i in range(n-1): # Heisenberg couplings
        A0 = A0 + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(n): # terms that make A non-Hermitian
        A0 = A0 + 0.3j*sc.Sy[i]*sc.Sx[i]
    A = sc.toMPO(A0) # an operator that can be multiplied
    Ad = A.get_dagger()
    tr["A"].append(abs(A.trace()))
    tr["A*A"].append(abs((A*A).trace()))
    tr["A*A^dagger"].append(abs((A*Ad).trace()))
    tr["A*A^dagger*A"].append(abs((A*Ad*A).trace()))
    print("n =", n, "Tr(A*A^dagger) =", (A*Ad).trace())

import matplotlib.pyplot as plt
for key in tr:
    plt.plot(ns, tr[key], marker="o", label=key)
plt.xlabel("chain length n")
plt.ylabel("|Trace|")
plt.yscale("log")
plt.legend()
plt.show()
