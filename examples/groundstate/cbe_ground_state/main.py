# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Ground-state energy of the S=1/2 Heisenberg chain with two DMRG solvers:
# two-site DMRG, and one-site DMRG with controlled bond expansion (CBE).
# CBE adds to each bond the few states where H|psi> has weight outside
# the current bond space, so it converges like two-site DMRG while the
# local problem is a one-site one. We plot the error against ED.
import numpy as np
from dmrgpy import spinchain

n = 12 # number of sites
ms = [2, 4, 8, 16, 32] # bond dimensions
errors = {"TWO_SITE": [], "CBE": []}
for method in errors:
    for m in ms:
        np.random.seed(1) # same random start for both solvers
        sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)],
                                  itensor_version="python")
        h = 0
        for i in range(n-1):
            h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
        sc.set_hamiltonian(h)
        sc.maxm = m # bond dimension
        sc.nsweeps = 8 # number of sweeps
        sc.dmrg_method = method # "TWO_SITE" (default) or "CBE"
        errors[method].append(sc.gs_energy() - sc.gs_energy(mode="ED"))

print("Error at maxm=32, two-site:", errors["TWO_SITE"][-1])
print("Error at maxm=32, CBE:", errors["CBE"][-1])
assert abs(errors["CBE"][-1]) < 1e-6

import matplotlib.pyplot as plt
plt.semilogy(ms, np.abs(errors["TWO_SITE"]) + 1e-16, "o-", linewidth=4,
             label="two-site DMRG")
plt.semilogy(ms, np.abs(errors["CBE"]) + 1e-16, "x--", markersize=10,
             label="CBE")
plt.xlabel("bond dimension maxm")
plt.ylabel("E - E_ED")
plt.legend()
plt.show()
