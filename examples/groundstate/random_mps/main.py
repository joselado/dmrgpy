# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# A random matrix product state (MPS) is a normalized wavefunction,
# <psi|psi> = 1, whatever the size of the chain.
# It is used as a starting point of DMRG and of time evolution.
# We check the norm for several chain lengths.
from dmrgpy import spinchain

sizes = [4, 6, 8, 10, 14] # chain lengths
norms = []
for n in sizes:
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
    wf = sc.random_mps() # random wavefunction
    norms.append(wf.overlap(wf).real)
    print("n =", n, " <psi|psi> =", norms[-1])

import matplotlib.pyplot as plt
plt.plot(sizes, norms, "o-", label="random MPS")
plt.axhline(1.0, color="gray", linestyle="--", label="expected")
plt.xlabel("chain length n")
plt.ylabel("<psi|psi>")
plt.legend()
plt.show()
