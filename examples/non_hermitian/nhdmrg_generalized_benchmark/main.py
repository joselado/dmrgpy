# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Generalized eigenvalue problem H|psi> = lambda A|psi> for a non-Hermitian
# H, the interacting fermion chain with a staggered imaginary onsite term,
# and a Hermitian positive metric A = 1 + 0.5 sum_i N_i. The ground state
# is the eigenvalue of smallest real part. We solve it with the
# non-Hermitian DMRG and compare with the exact generalized eigenvalues.
import numpy as np
import scipy.linalg as sla
from dmrgpy import fermionchain

ns = [4, 6] # chain lengths
exact, dmrg = [], []
for n in ns:
    rng = np.random.RandomState(2)
    fc = fermionchain.Fermionic_Chain(n)
    fc.setup_python() # pure Python backend
    h = 0
    for i in range(n-1):
        h = h + rng.random()*(fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i])
        h = h + rng.random()*(fc.N[i]-0.5)*(fc.N[i+1]-0.5)
    for i in range(n):
        h = h + 1j*(-1)**i*0.3*fc.N[i] # non-Hermitian term
    a = 1 + 0.5*sum(fc.N) # metric
    fc.set_hamiltonian(h)
    fc.maxm = 60 # bond dimension
    fc.nsweeps = 20 # number of DMRG sweeps
    dmrg.append(fc.gs_energy_generalized(a))
    # exact solution from the matrices
    ed = fc.get_ED_obj()
    hm = ed.MO2matrix(h).toarray()
    am = ed.MO2matrix(a).toarray()
    w = sla.eig(hm, am, right=False)
    exact.append(w[np.argsort(w.real)][0])
print("Max difference =", np.max(np.abs(np.array(exact)-dmrg)))

import matplotlib.pyplot as plt
plt.plot(ns, np.real(exact), "o-", linewidth=4, label="exact")
plt.plot(ns, np.real(dmrg), "x--", linewidth=2, label="NH-DMRG")
plt.xlabel("number of sites n")
plt.ylabel("Re(lambda)")
plt.legend()
plt.show()
