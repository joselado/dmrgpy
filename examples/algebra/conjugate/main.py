# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain

def get_dots(n):
    """Build a Heisenberg chain of length n, compute <wf0|wf0> (should be
    1) and <wf0|i*wf0*> (should be -i), and return both as complex
    numbers."""
    spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
    sc = spinchain.Spin_Chain(spins) # create the spin chain
    h = 0 # generate a Heisenberg Hamiltonian
    for i in range(n-1):
        h = h +sc.Sx[i]*sc.Sx[i+1]
        h = h +sc.Sy[i]*sc.Sy[i+1]
        h = h +sc.Sz[i]*sc.Sz[i+1]

    # the dagger of an operator con be computed as
    # B = A.get_dagger()

    sc.set_hamiltonian(h) # set the Hamiltonian
    wf0 = sc.get_gs(mode="DMRG") # compute ground state
    wf0i = 1j*wf0 # imaginary unit times ground state

    wf1 = wf0i.get_conjugate() # compute the conjugate

    return wf0.dot(wf0), wf0.dot(wf1)

ns = [2,4,6,8] # sweep the chain length
norms,confs = [],[]
for n in ns:
    print("Chain length",n)
    norm,conf = get_dots(n)
    print("This should be 1")
    print(norm)
    print("This should be -i")
    print(conf)
    norms.append(norm)
    confs.append(conf)

import matplotlib.pyplot as plt
fig,axes = plt.subplots(1,2,figsize=(10,4))
axes[0].plot(ns,[abs(x) for x in norms],marker="o",label="$|\\langle \\psi_0|\\psi_0\\rangle|$")
axes[0].axhline(1.0,color="gray",linestyle="--",label="expected (1)")
axes[0].set_xlabel("Chain length n")
axes[0].set_ylabel("$|\\langle \\psi_0|\\psi_0\\rangle|$")
axes[0].legend()

axes[1].plot(ns,[x.real for x in confs],marker="o",label="Re")
axes[1].plot(ns,[x.imag for x in confs],marker="s",label="Im")
axes[1].axhline(-1.0,color="gray",linestyle="--",label="expected Im (-1)")
axes[1].set_xlabel("Chain length n")
axes[1].set_ylabel("$\\langle \\psi_0|\\psi_1\\rangle$ (conjugate overlap)")
axes[1].legend()

plt.tight_layout()
plt.show()
