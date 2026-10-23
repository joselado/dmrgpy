# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Generalized eigenvalue problem H|psi> = lambda A|psi>, with H an
# interacting fermion chain with random hoppings and A = 1 + 0.5 sum_i n_i
# a positive metric. It appears when the basis is not orthogonal.
# We solve it with DMRG (gs_energy_generalized) and compare with the
# exact solution of scipy for several chain lengths.
import numpy as np
import scipy.linalg as sla
from dmrgpy import fermionchain

lam, exact, ns = [], [], [6, 7, 8]
for n in ns:
    np.random.seed(2)
    fc = fermionchain.Fermionic_Chain(n, itensor_version="python")
    h = 0
    for i in range(n-1):
        h = h + np.random.random()*(fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i])
        h = h + np.random.random()*(fc.N[i]-.5)*(fc.N[i+1]-.5)
    a = 1 # metric
    for i in range(n): a = a + 0.5*fc.N[i]
    fc.set_hamiltonian(h)
    fc.maxm = 60 # bond dimension
    fc.nsweeps = 14 # number of DMRG sweeps
    fc.cutoff = 1e-12 # no truncation
    lam.append(fc.gs_energy_generalized(a))
    ed = fc.get_ED_obj() # exact solution
    Hm, Am = ed.MO2matrix(h).toarray(), ed.MO2matrix(a).toarray()
    exact.append(sla.eigh(Hm, Am, eigvals_only=True)[0])
print("Max difference =", np.max(np.abs(np.array(lam)-exact)))

import matplotlib.pyplot as plt
plt.plot(ns, exact, linewidth=4, label="exact")
plt.plot(ns, lam, "o--", label="DMRG")
plt.xlabel("number of sites n")
plt.ylabel("lowest eigenvalue")
plt.legend()
plt.show()
