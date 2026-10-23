# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Mutual information I(A:B) = S_A + S_B - S_AB between sites 1 and 2
# and sites 3 and 4 of a Heisenberg chain with four spins. We weaken
# the coupling in the middle and see how the correlation between the
# two halves disappears. It measures the total correlation of two parts.
import numpy as np
from dmrgpy import spinchain

n = 4 # number of sites

def get(c01):
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
    h = 0
    for i in range(n-1):
        c = 1.0
        if i == 1: c = c01 # coupling in the middle
        h = h + c*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
    sc.set_hamiltonian(h)
    wf = sc.get_gs() # ground state
    return wf.get_mutual_information(1, 2)

cs = np.linspace(0., 1., 10) # coupling between the two halves
ss = [get(c) for c in cs]

import matplotlib.pyplot as plt
plt.plot(cs, ss, marker="o")
plt.xlabel("coupling between (1,2) and (3,4)")
plt.ylabel("mutual information")
plt.show()
