# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Complex conjugate of a many-body wavefunction. For the ground state
# psi0 of a real Hamiltonian, <psi0|psi0> = 1. If we take psi1 = conj(i*psi0),
# which is -i*psi0, then <psi0|psi1> = -i. This is a check that the
# conjugation acts on the coefficients of the state.
import numpy as np
from dmrgpy import spinchain

ns = [2,4,6,8] # chain lengths
norms, overlaps = [], []
for n in ns:
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
    h = 0
    for i in range(n-1): # Heisenberg Hamiltonian
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = 20 # bond dimension
    sc.nsweeps = 10 # number of DMRG sweeps
    wf0 = sc.get_gs(mode="DMRG") # ground state
    wf1 = (1j*wf0).get_conjugate() # conjugate of i*psi0
    norms.append(wf0.dot(wf0))
    overlaps.append(wf0.dot(wf1))
print("Overlap <psi0|psi1> (should be -i):", overlaps)

import matplotlib.pyplot as plt
plt.plot(ns, np.real(overlaps), "o-", label="Re <psi0|psi1>")
plt.plot(ns, np.imag(overlaps), "s-", label="Im <psi0|psi1>")
plt.plot(ns, np.real(norms), "^-", label="<psi0|psi0>")
plt.axhline(-1, color="gray", linestyle="--")
plt.xlabel("chain length n")
plt.ylabel("overlap")
plt.legend()
plt.show()
