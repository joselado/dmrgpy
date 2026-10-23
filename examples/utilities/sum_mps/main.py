# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Algebra with matrix product states (MPS) on a chain of 20 spins. The
# ground states of +-sum Sz are the fully up and fully down states, which
# are orthogonal. Applying 2 Sx on every site flips all the spins, and the
# sum of two MPS is again an MPS. We follow the overlaps as a function of n.
import numpy as np
from dmrgpy import spinchain

def overlaps(n):
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
    hfe = sum(sc.Sz) # fully ferromagnetic operator
    sc.set_hamiltonian(-hfe) ; wfup = sc.get_gs() # all spins up
    sc.set_hamiltonian(hfe) ; wfdn = sc.get_gs() # all spins down
    wfflip = wfdn.copy()
    for i in range(n):
        wfflip = 2.*sc.Sx[i]*wfflip # flip every spin
    wfs = wfup + 3.7*wfdn # sum of two MPS
    return wfup.dot(wfdn).real, wfup.dot(wfflip).real, wfdn.overlap(wfs).real

ns = range(4, 21, 4) # chain sizes
out = np.array([overlaps(n) for n in ns])
print("Overlaps for the largest chain =", out[-1])

import matplotlib.pyplot as plt
plt.plot(ns, out[:,0], "o-", label="up and down")
plt.plot(ns, out[:,1], "o-", label="up and flipped down")
plt.plot(ns, out[:,2], "o-", label="down and (up + 3.7 down)")
plt.xlabel("chain size n")
plt.ylabel("overlap")
plt.legend()
plt.show()
