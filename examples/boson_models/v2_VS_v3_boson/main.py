# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Ground state energy of a Bose-Hubbard chain with random hoppings,
# H = sum_i t_i (b_i^dagger b_{i+1} + h.c.) + U sum_i n_i^2.
# We compare the ITensor v2, v3 and pure Python DMRG backends with
# exact diagonalization (ED), as a function of the interaction U.
import numpy as np
from dmrgpy import bosonchain

n = 6 # number of sites

def energy(U, version, mode="DMRG"):
    bc = bosonchain.Bosonic_Chain(n)
    if version == "python": bc.setup_python()
    else: bc.setup_cpp(version)
    bc.maxm = 100 # bond dimension, no truncation here
    bc.nsweeps = 80 # number of DMRG sweeps
    np.random.seed(11) # same hoppings in all backends
    h = 0
    for i in range(n-1):
        h = h + np.random.random()*(bc.Adag[i]*bc.A[i+1] + bc.Adag[i+1]*bc.A[i])
    for i in range(n):
        h = h + U*bc.Adag[i]*bc.A[i]*bc.Adag[i]*bc.A[i]
    bc.set_hamiltonian(h)
    return bc.gs_energy(mode=mode)

Us = [0.1, 0.3, 0.5]
e2 = [energy(U, 2) for U in Us]
e3 = [energy(U, 3) for U in Us]
ep = [energy(U, "python") for U in Us]
ed = [energy(U, 2, mode="ED") for U in Us]
err = np.max(np.abs(np.array(e2+e3+ep)-np.array(ed*3)))
print("Max difference =", err)
tol = 1e-5 # tolerance of the check
assert err < tol

import matplotlib.pyplot as plt
plt.plot(Us, ed, linewidth=4, label="ED")
plt.plot(Us, e2, "o--", label="v2")
plt.plot(Us, e3, "x--", label="v3")
plt.plot(Us, ep, "+--", label="python")
plt.xlabel("U")
plt.ylabel("ground state energy")
plt.legend()
plt.show()
