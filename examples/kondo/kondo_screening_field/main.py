# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Kondo screening of an Anderson impurity in a magnetic field, with the
# NRG-MPS method of arXiv:2501.10746. The band is mapped to a Wilson chain
# and DMRG gives the ground state. We project it on a singly occupied
# impurity and compute its two Schmidt coefficients: both are 1/sqrt(2) for
# the Kondo singlet, and the larger one goes to 1 when a field B above the
# Kondo temperature polarizes the impurity.
import numpy as np
from dmrgpy import wilsonchain

U, tD, Lambda = 1.0, 0.3, 3.0 # repulsion, hybridization, discretization
n = 13 # Wilson chain sites, odd so that the ground state is a singlet
TK = np.sqrt(np.pi*tD**2/4)*np.exp(-U/(4*tD**2)) # estimate, flat band, D=U

Bs = np.logspace(-4, -0.5, 8) # magnetic fields
l1, l2 = [], []
for B in Bs:
    fc = wilsonchain.anderson_chain(n, ed=-U/2, U=U, tD=tD, Lambda=Lambda, B=B)
    fc.maxm = 100 # bond dimension
    fc.nsweeps = 12 # number of DMRG sweeps
    fc.cutoff = 1e-12
    r = fc.get_kondo_screening()
    l1.append(r["schmidt"][0])
    l2.append(r["schmidt"][1])
    print("B =", B, "Schmidt coefficients", r["schmidt"])

import matplotlib.pyplot as plt
plt.semilogx(Bs, l1, "o-", label="lambda_1")
plt.semilogx(Bs, l2, "o-", label="lambda_2")
plt.axhline(1/np.sqrt(2), color="gray", linewidth=0.5)
plt.axvline(TK, color="gray", linestyle="--", label="T_K estimate")
plt.xlabel("magnetic field B")
plt.ylabel("Schmidt coefficients")
plt.legend()
plt.show()
