# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import spinchain

# Sweep the chain length and check that a fresh random MPS is always
# normalized (overlap(wf,wf)==1) regardless of system size, and that the
# spin operator algebra (commutation relations) passes sc.test() at every
# size -- a cheap sweep instead of a single-size, single-scalar check.
sizes = [4, 6, 8, 10, 14]
overlaps = []
for n in sizes:
    spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
    sc = spinchain.Spin_Chain(spins) # create the spin chain
    wf = sc.random_mps() # random wavefunction
    ov = wf.overlap(wf).real
    print("n =",n,"overlap(wf,wf) =",ov)
    overlaps.append(ov)
    sc.test()

plt.plot(sizes, overlaps, "o-")
plt.axhline(1.0, color="k", linestyle="--", label="expected norm")
plt.xlabel("chain length n")
plt.ylabel(r"$\langle\psi|\psi\rangle$ of a random MPS")
plt.title("Random MPS normalization vs system size")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()











